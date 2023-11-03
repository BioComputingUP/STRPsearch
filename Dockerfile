# Use a base image with Conda and Python
FROM continuumio/miniconda3

# Set the working directory inside the container
WORKDIR /app

# Create a conda environment and install dependencies
RUN conda create --name bio38_env python=3.8

## Activate the environment
RUN echo "source activate bio38_env" > ~/.bashrc
ENV PATH /opt/conda/envs/bio38_env/bin:$PATH

# Install the required packages using conda
RUN conda install -c bioconda tmalign
RUN conda install -c conda-forge -c bioconda foldseek

# Install Python dependencies
RUN pip install biopython numpy pandas matplotlib requests scipy

# Copy the entire project directory into the container
COPY . /app

# Define the command to run your script
ENTRYPOINT ["python3", "bin/repeatsdb-lite-2.py"]
