# Step 1: Base image
FROM ubuntu:20.04

# Step 2: Set noninteractive mode (prevents UI prompts)
ENV DEBIAN_FRONTEND=noninteractive

# Step 3: Install basic tools and dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    libgsl-dev \
    zlib1g-dev \
    pkg-config \
    bzip2 \
    ca-certificates \
    libglib2.0-0 \
    libxext6 \
    libsm6 \
    libxrender1 \
    gcc \
    && rm -rf /var/lib/apt/lists/*

# Step 4: Install Anaconda
ENV CONDA_DIR=/opt/conda
ENV PATH=$CONDA_DIR/bin:$PATH

RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2023.07-1-Linux-x86_64.sh -O /tmp/anaconda.sh && \
    bash /tmp/anaconda.sh -b -p $CONDA_DIR && \
    rm /tmp/anaconda.sh && \
    conda clean -afy
RUN conda init bash

# Step 5: Set work directory
WORKDIR /app 

# Step 6: Copy your scripts (optional if needed)
# COPY . /app 
COPY . .
#RUN conda env create --file=/app/env/ScanFold3.yml
#RUN conda activate ScanFold3
#RUN pip install --no-cache-dir -r requirements.txt

# Step 7: Default command (can be overridden)
CMD [ "bash" ]

 
