FROM python:3.11

RUN pip install protmapper[rest_api] && \
    python -m protmapper.resources

ENTRYPOINT python -m protmapper.rest_api
