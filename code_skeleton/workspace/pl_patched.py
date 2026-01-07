import os
import pickle
import lmdb
import torch
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from torch.utils.data import Dataset
from tqdm.auto import tqdm

from ..protein_ligand import PDBProtein, parse_sdf_file
from ..data import ProteinLigandData, torchify_dict

def process_one_item(args):
    """
    Worker function to process a single data item.
    Returns (key_bytes, pickled_value_bytes) or None if skipped.
    """
    i, (pocket_fn, ligand_fn, _, rmsd_str), raw_path = args
    if pocket_fn is None:
        return None
        
    try:
        pocket_dict = PDBProtein(os.path.join(raw_path, pocket_fn)).to_dict_atom()
        ligand_dict = parse_sdf_file(os.path.join(raw_path, ligand_fn))
        data = ProteinLigandData.from_protein_ligand_dicts(
            protein_dict=torchify_dict(pocket_dict),
            ligand_dict=torchify_dict(ligand_dict),
        )
        data.protein_filename = pocket_fn
        data.ligand_filename = ligand_fn
        
        # Determine ID based on index 'i' passed in
        key = str(i).encode()
        value = pickle.dumps(data)
        return (key, value)
    except Exception as e:
        # In multiprocessing, simple print might interleave, but it's okay for debugging
        print(f'Skipping ({i}) {ligand_fn} for: {e}')
        return None

class PocketLigandPairDataset(Dataset):

    def __init__(self, raw_path, transform=None):
        super().__init__()
        self.raw_path = raw_path.rstrip('/')
        self.index_path = os.path.join(self.raw_path, 'index.pkl')
        self.processed_path = os.path.join(os.path.dirname(self.raw_path), os.path.basename(self.raw_path) + '_processed.lmdb')
        self.name2id_path = os.path.join(os.path.dirname(self.raw_path), os.path.basename(self.raw_path) + '_name2id.pt')
        self.transform = transform
        self.db = None

        self.keys = None

        if not (os.path.exists(self.processed_path) and os.path.exists(self.name2id_path)):
            self._process()
            self._precompute_name2id()

        self.name2id = torch.load(self.name2id_path)

    def _connect_db(self):
        """
            Establish read-only database connection
        """
        assert self.db is None, 'A connection has already been opened.'
        self.db = lmdb.open(
            self.processed_path,
            map_size=10*(1024*1024*1024),   # 10GB
            create=False,
            subdir=False,
            readonly=True,
            lock=False,
            readahead=False,
            meminit=False,
        )
        with self.db.begin() as txn:
            self.keys = list(txn.cursor().iternext(values=False))

    def _close_db(self):
        self.db.close()
        self.db = None
        self.keys = None
        
    def _process(self):
        print(f"Creating LMDB database at {self.processed_path}...")
        db = lmdb.open(
            self.processed_path,
            map_size=10*(1024*1024*1024),   # 10GB
            create=True,
            subdir=False,
            readonly=False, # Writable
        )
        with open(self.index_path, 'rb') as f:
            index = pickle.load(f)

        # Multiprocessing Setup
        num_workers = max(1, multiprocessing.cpu_count() - 1)
        print(f"Using {num_workers} workers for data processing.")
        
        tasks = []
        for i, item in enumerate(index):
            tasks.append((i, item, self.raw_path))

        # Process in parallel
        # We hold results in memory to write in one transaction, or we could write in batches.
        # Given dataset size (~10k items), holding in memory is safe and simpler.
        
        valid_count = 0
        with db.begin(write=True, buffers=True) as txn:
             with ProcessPoolExecutor(max_workers=num_workers) as executor:
                # Use executor.map to preserve order (though for LMDB keys it matters less if we sort, but 'i' is key)
                results = list(tqdm(executor.map(process_one_item, tasks), total=len(tasks), desc="Processing SDFs"))
                
                for res in results:
                    if res is not None:
                        key, val = res
                        txn.put(key, val)
                        valid_count += 1
                        
        print(f"Written {valid_count} items to LMDB.")
        db.close()

    def _precompute_name2id(self):
        name2id = {}
        for i in tqdm(range(self.__len__()), 'Indexing'):
            try:
                data = self.__getitem__(i)
            except AssertionError as e:
                print(i, e)
                continue
            name = (data.protein_filename, data.ligand_filename)
            name2id[name] = i
        torch.save(name2id, self.name2id_path)
    
    def __len__(self):
        if self.db is None:
            self._connect_db()
        return len(self.keys)

    def __getitem__(self, idx):
        if self.db is None:
            self._connect_db()
        key = self.keys[idx]
        data = pickle.loads(self.db.begin().get(key))
        data.id = idx
        assert data.protein_pos.size(0) > 0
        if self.transform is not None:
            data = self.transform(data)
        return data
        

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('path', type=str)
    args = parser.parse_args()

    PocketLigandPairDataset(args.path)
