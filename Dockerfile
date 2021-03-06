FROM python:2.7.17-slim AS base

FROM base AS builder

# Build dependencies
RUN apt-get update -y
RUN apt-get install -y git
RUN apt-get install -y pkg-config build-essential liblzma-dev

# Python dependencies
COPY requirements.txt /tmp/requirements.txt
RUN pip install -r /tmp/requirements.txt
# Optional dependency
RUN pip install pyliblzma

# The repository itself
RUN mkdir -p /tmp/agora
COPY . /tmp/agora
WORKDIR /tmp/agora
RUN git clean -d -f
RUN git clean -X -f
RUN rm -rf .git

# Compile all Python files to speed things up
RUN python -m compileall .

# Trick to only keep the files we need and make the image smaller
FROM base AS prod
COPY --from=builder /tmp/agora /opt/agora
COPY --from=builder /usr/local /usr/local

# Test the image
FROM prod AS test
WORKDIR /opt/agora
RUN ./checkAgoraIntegrity.sh

FROM prod
