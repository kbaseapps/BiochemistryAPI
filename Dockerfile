FROM kbase/kbase:sdkbase.latest
MAINTAINER KBase Developer
# -----------------------------------------

# Here we install a python coverage tool and an
# https library that is out of date in the base image.

# update security libraries in the base image
RUN pip install setuptools --upgrade \
    && pip install cffi --upgrade \
    && pip install pyopenssl --upgrade \
    && pip install ndg-httpsclient --upgrade \
    && pip install pyasn1 --upgrade \
    && pip install requests --upgrade \
    && pip install 'requests[security]' --upgrade

# -----------------------------------------
# download and use conda to install rdkit and add the installed libs to the paths
RUN echo '$PATH:export PATH=/opt/conda/bin' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.3.27-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh
ENV PATH $PATH:/opt/conda/bin
RUN conda config --add channels  https://conda.anaconda.org/rdkit && \
    conda install -y nose \
                     cairo \
                     nomkl \
                     rdkit
ENV PYTHONPATH $PYTHONPATH:/usr/lib/python2.7/:/opt/conda/lib/python2.7/site-packages/
ENV LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu/:/opt/conda/lib/

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
