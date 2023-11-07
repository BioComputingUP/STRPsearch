# Use a base image with Conda and Python
FROM continuumio/miniconda3

# Set the working directory inside the container
WORKDIR /app

# Copy the entire project directory into the container
COPY . /app

# Create a conda environment and install dependencies
RUN conda env create -f rdblite_env.yml

## Activate the environment
RUN echo "source activate rdblite_env" > ~/.bashrc
ENV PATH /opt/conda/envs/rdblite_env/bin:$PATH

# Define the command to run your script
ENTRYPOINT ["python3", "bin/main.py"]
