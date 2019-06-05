FROM python:3.6

RUN pip install protmapper && \
    python -m protmapper.resources

ENTRYPOINT python -m protmapper.rest_api
