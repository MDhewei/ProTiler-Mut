# Use a slim Python base
FROM continuumio/miniconda3

# Set working directory
WORKDIR /app

# Install system dependencies (if needed)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Install pymol 
RUN conda install -c conda-forge pymol-open-source

# Copy your source code and project files
COPY . .

# Install your package
RUN pip install --upgrade pip && \
    pip install .

# Set CLI entrypoint
ENTRYPOINT ["protiler-mut"]
CMD ["--help"]
