# Use the official pre-built gnina image from Docker Hub.
FROM gnina/gnina:latest

# Set the working directory
WORKDIR /app

# Copy the docking script into the container directory /usr/local/bin
COPY run_gnina_docking.sh /usr/local/bin/
COPY run_gnina_redock.sh /usr/local/bin/

# Make the script executable
RUN chmod +x /usr/local/bin/run_gnina_docking.sh /usr/local/bin/run_gnina_redock.sh
