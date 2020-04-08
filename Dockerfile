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

FROM base AS prod
ARG INSTALLATION_PATH=/opt/agora
COPY --from=builder /tmp/agora "$INSTALLATION_PATH"
COPY --from=builder /usr/local /usr/local

FROM prod AS test
WORKDIR /opt/agora
RUN ./checkAgoraIntegrity.sh

FROM prod
