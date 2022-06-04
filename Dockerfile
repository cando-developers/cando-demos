FROM ubuntu:jammy

SHELL ["/bin/bash", "-c"]

ARG NB_USER=app
ARG NB_UID=1000

ENV DEBIAN_FRONTEND=noninteractive
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}
ENV PATH "${HOME}/.local/bin:${PATH}"

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -o Dpkg::Options::="--force-overwrite" -y \
        nano wget sudo git locales curl python3-pip nodejs npm && \
    bash -c "$(curl -fsSL https://www.thirdlaw.tech/pkg/cando.sh)"

RUN echo 'en_US.UTF-8 UTF-8' >/etc/locale.gen
RUN sudo -E locale-gen

RUN useradd --create-home --shell=/bin/false --uid=${NB_UID} ${NB_USER} && \
    usermod -aG sudo $NB_USER && \
    passwd -d $NB_USER

WORKDIR ${HOME}

COPY . demos/
RUN chown -R ${NB_UID} demos/ && \
    chgrp -R ${NB_USER} demos/

WORKDIR ${HOME}/demos/
USER ${NB_USER}

RUN pip install notebook jupyterlab && \
    jupyter-labextension install cytoscape-clj kekule-clj ngl-clj \
      resizable-box-clj @jupyter-widgets/jupyterlab-manager \
      jupyterlab_templates jupyterlab-debugger-restarts jupyterlab-molviewer
