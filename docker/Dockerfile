ARG FROM_IMAGE=gitlab-registry.cern.ch/cms-cloud/cmssw-docker/al9-cms:latest
FROM --platform=linux/amd64 ${FROM_IMAGE}

SHELL ["/bin/bash", "-c"]

ARG CMSSW_VERSION
ARG COMBINE_TAG
ARG HARVESTER_BRANCH

USER cmsusr
WORKDIR /home/cmsusr

ADD --chown=cmsusr:cmsusr ${CMSSW_VERSION} /home/cmsusr/${CMSSW_VERSION}
RUN mkdir -p "/home/cmsusr/${CMSSW_VERSION}/src/2DAlphabet"

COPY --chown=cmsusr:cmsusr TwoDAlphabet "/home/cmsusr/${CMSSW_VERSION}/src/2DAlphabet/TwoDAlphabet"
COPY --chown=cmsusr:cmsusr test "/home/cmsusr/${CMSSW_VERSION}/src/2DAlphabet/test"
COPY --chown=cmsusr:cmsusr setup.py example.py example_config.json "/home/cmsusr/${CMSSW_VERSION}/src/2DAlphabet/"
COPY --chown=cmsusr:cmsusr docker/setup.sh /home/cmsusr/setup_2dalphabet.sh

ENV CMSSW_VERSION=${CMSSW_VERSION}
ENV COMBINE_TAG=${COMBINE_TAG}
ENV HARVESTER_BRANCH=${HARVESTER_BRANCH}
ENV TWODALPHABET_PATH=/home/cmsusr/${CMSSW_VERSION}/src/2DAlphabet
ENV PATH=/home/cmsusr/.local/bin:${PATH}

RUN source /cvmfs/cms.cern.ch/cmsset_default.sh && \
    cd "/home/cmsusr/${CMSSW_VERSION}/src" && \
    scram b ProjectRename && \
    eval "$(scramv1 runtime -sh)" && \
    scram b -j"$(nproc)" && \
    python3 -m pip install --user --no-cache-dir virtualenv && \
    python3 -m virtualenv /home/cmsusr/twoD-env && \
    source /home/cmsusr/twoD-env/bin/activate && \
    python3 -m pip install --no-cache-dir -e "${TWODALPHABET_PATH}" && \
    python3 -c "from ROOT import RooParametricHist2D; import TwoDAlphabet"

RUN chmod +x /home/cmsusr/setup_2dalphabet.sh && \
    echo "source /home/cmsusr/setup_2dalphabet.sh" >> /home/cmsusr/.bashrc

CMD ["/bin/bash"]
