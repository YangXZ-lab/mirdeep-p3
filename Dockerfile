# mirdeep-p3 Docker image
# Build:  docker build -t mirdeep-p3:3.1.4c .
# Run:    docker run --rm -v $(pwd):/data mirdeep-p3:3.1.4c mirdeep-p3 -h

FROM continuumio/miniconda3:latest

# ------------------------------------------------------------------
# 1. System tools needed by conda + pipeline
# ------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    wget \
    curl \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------------
# 2. Create conda env from the project's environment file
#    (uses official conda-forge/bioconda channels, not mirrors)
# ------------------------------------------------------------------
COPY mirdp3_environment.yml /tmp/mirdp3_environment.yml
RUN conda env create -f /tmp/mirdp3_environment.yml -n mirdp3 \
    && conda clean -afy

# ------------------------------------------------------------------
# 3. Make the env the default when running the container
# ------------------------------------------------------------------
ENV PATH=/opt/conda/envs/mirdp3/bin:$PATH
RUN echo "source activate mirdp3" > /root/.bashrc

# ------------------------------------------------------------------
# 4. Copy software source tree (mirdeep-p3 locates src/ relative to itself)
# ------------------------------------------------------------------
WORKDIR /opt/mirdeep-p3
COPY . /opt/mirdeep-p3/

RUN chmod 755 /opt/mirdeep-p3/mirdeep-p3

# ------------------------------------------------------------------
# 5. Default workdir for user data (mount with -v)
# ------------------------------------------------------------------
WORKDIR /data

# Default entry: print help if no args given
ENTRYPOINT ["/opt/mirdeep-p3/mirdeep-p3"]
CMD ["-h"]
