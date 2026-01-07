FROM gnina/gnina:latest

# Screening Script (Raw SMILES input)
COPY run_gnina_docking.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/run_gnina_docking.sh

ENTRYPOINT ["/bin/bash"]
