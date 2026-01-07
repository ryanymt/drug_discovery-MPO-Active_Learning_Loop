FROM gnina/gnina:latest

# Redocking Script (CSV input + SDF output)
COPY run_gnina_redock.sh /usr/local/bin/
RUN chmod +x /usr/local/bin/run_gnina_redock.sh

ENTRYPOINT ["/bin/bash"]
