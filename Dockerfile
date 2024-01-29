# Use a base image with Conda and Python
FROM continuumio/miniconda3

# Set the working directory inside the container
WORKDIR /app

# Copy the entire project directory into the container
COPY . /app

# Create a conda environment and install dependencies
RUN conda env create -f environment.yml

## Activate the environment
RUN echo "source activate strpsearch_env" > ~/.bashrc
ENV PATH /opt/conda/envs/strpsearch_env/bin:$PATH

RUN apt-get update && apt-get install ffmpeg libsm6 libxext6  -y

# Define the command to run your script
ENTRYPOINT ["python3", "bin/strpsearch.py"]
