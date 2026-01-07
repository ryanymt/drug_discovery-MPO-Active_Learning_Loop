import argparse
import os
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM
import pandas as pd
from huggingface_hub import login

def get_prompt(smiles):
    """Creates the prompt for the TxGemma model using TDC format for Toxicity Classification."""
    return (
        f"<start_of_turn>user\n"
        f"Instructions: Answer the following question about drug properties.\n"
        f"Context: Clinical toxicity is a major cause of attrition in drug development. We need to identify if a molecule is toxic.\n"
        f"Question: Given a drug SMILES string, predict whether it is (A) Non-toxic or (B) Toxic. Drug SMILES: {smiles}\n"
        f"Answer:<end_of_turn>\n"
        f"<start_of_turn>model\n"
    )

def run_inference(input_file, output_file):
    """
    Runs TxGemma inference on a file of SMILES strings and saves the predictions.
    """
    # --- Hugging Face Authentication ---
    hf_token = os.getenv("HF_TOKEN")
    if not hf_token:
        raise ValueError("Hugging Face token not found in environment variable HF_TOKEN")
    login(token=hf_token)
    print("Successfully logged in to Hugging Face Hub.")

    # --- Model and Tokenizer Loading ---
    model_id = "google/txgemma-9b-predict"
    
    if not torch.cuda.is_available():
        raise RuntimeError("This model requires a GPU, but CUDA is not available.")
    
    print(f"Loading model: {model_id} with 4-bit quantization.")
    tokenizer = AutoTokenizer.from_pretrained(model_id)
    model = AutoModelForCausalLM.from_pretrained(
        model_id,
        torch_dtype=torch.bfloat16,
        device_map="auto",
        trust_remote_code=True,
        load_in_4bit=True # Enable 4-bit quantization
    )
    print("Model and tokenizer loaded successfully.")

    # --- Data Loading ---
    try:
        with open(input_file, 'r') as f:
            smiles_list = [line.strip() for line in f if line.strip()]
    except FileNotFoundError:
        print(f"Error: Input file not found at {input_file}")
        return
        
    if not smiles_list:
        print("Input file is empty. Exiting.")
        pd.DataFrame().to_csv(output_file, index=False)
        return

    print(f"Processing {len(smiles_list)} SMILES strings from {input_file}")

    # --- Inference ---
    print("Generating predictions...")
    
    batch_size = 8
    all_results = []
    
    total_smiles = len(smiles_list)
    
    for i in range(0, total_smiles, batch_size):
        batch_smiles = smiles_list[i : i + batch_size]
        print(f"Processing batch {i // batch_size + 1}/{(total_smiles + batch_size - 1) // batch_size} (Indices {i} to {min(i + batch_size, total_smiles)})")
        
        prompts = [get_prompt(s) for s in batch_smiles]
        inputs = tokenizer(prompts, return_tensors="pt", padding=True).to(model.device)
        
        input_length = inputs.input_ids.shape[1]
        
        with torch.no_grad():
            outputs = model.generate(**inputs, max_new_tokens=256, do_sample=False)
        
        # Slice to keep only the generated tokens
        generated_tokens = outputs[:, input_length:]
        
        # Decode only the generation
        decoded_outputs = tokenizer.batch_decode(generated_tokens, skip_special_tokens=True)
        
        # --- Results Parsing for Batch ---
        for j, prediction_text in enumerate(decoded_outputs):
            prediction_text = prediction_text.strip()
            print(f"RAW OUTPUT for SMILES {batch_smiles[j]}: {prediction_text!r}")
            
            properties = {}
            properties['raw_prediction'] = prediction_text
            
            # Parse (A)/(B) output
            if '(A)' in prediction_text or '(B)' in prediction_text:
                if '(A)' in prediction_text:
                    properties['toxicity_pred'] = 0  # Non-toxic
                    properties['toxicity_label'] = 'Non-toxic'
                elif '(B)' in prediction_text:
                    properties['toxicity_pred'] = 1  # Toxic
                    properties['toxicity_label'] = 'Toxic'
            else:
                 # Fallback: check for keywords if A/B is missing
                lower_text = prediction_text.lower()
                if 'non-toxic' in lower_text:
                     properties['toxicity_pred'] = 0
                     properties['toxicity_label'] = 'Non-toxic'
                elif 'toxic' in lower_text:
                     properties['toxicity_pred'] = 1
                     properties['toxicity_label'] = 'Toxic'
                else:
                    properties['toxicity_pred'] = None
                    properties['toxicity_label'] = 'Unknown'

            properties['smiles'] = batch_smiles[j]
            all_results.append(properties)

    df = pd.DataFrame(all_results)
    
    # Rename 'toxicity_pred' to 'toxicity' for cleanliness
    if 'toxicity_pred' in df.columns:
        df = df.rename(columns={'toxicity_pred': 'toxicity'})
    
    # Final output columns: Just smiles and toxicity (0/1)
    final_cols = ['smiles', 'toxicity']
    
    # Filter to ensure columns exist (fallback if something went wrong, though init guarantees it)
    # If 'toxicity' is missing (e.g. no predictions parsed?), fill with None or -1? 
    # For now, just save what we have that matches.
    df = df[[c for c in final_cols if c in df.columns]]
    
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Run TxGemma ADMET prediction.")
    parser.add_argument("--input-file", required=True, help="Path to the input file with SMILES strings.")
    parser.add_argument("--output-file", required=True, help="Path to save the output CSV file.")
    
    args = parser.parse_args()
    
    run_inference(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
