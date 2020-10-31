FROM thirdlaw/cando:latest

ARG NB_USER=app
ARG NB_UID=1000

ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

RUN rm -rf /opt/clasp/run/kernels/cando/
COPY kernels /opt/clasp/run/kernels/
COPY . ${HOME}/demos/
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
WORKDIR ${HOME}/demos/

ENTRYPOINT ["scripts/docker-entry.sh"]
