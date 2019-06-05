FROM python:3.6

RUN pip install git+https://github.com/bgyori/protmapper.git@docker#egg=protmapper[rest_api] && \
    python -m protmapper.resources

ENTRYPOINT python -m protmapper.rest_api
