FROM python:2.7.17-slim AS base

FROM base AS builder

# Build dependencies
RUN apt-get update -y
RUN apt-get install -y git

# The repository itself
RUN mkdir -p /tmp/agora
COPY . /tmp/agora
WORKDIR /tmp/agora
RUN git clean -d -f
RUN git clean -X -f
RUN rm -rf .git

# Python dependencies
RUN pip install -r requirements.txt

FROM base AS prod
ARG INSTALLATION_PATH=/opt/agora
COPY --from=builder /tmp/agora "$INSTALLATION_PATH"
COPY --from=builder /usr/local /usr/local
ENV PYTHONPATH=$INSTALLATION_PATH/LibsDyoGen/

FROM prod AS test
WORKDIR /opt/agora
#RUN ./checkAgoraIntegrity.sh

FROM prod
