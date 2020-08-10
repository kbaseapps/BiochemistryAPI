FROM kbase/sdkbase2:python
MAINTAINER KBase Developer2

# -----------------------------------------
RUN apt-get update
RUN conda config --add channels  https://conda.anaconda.org/rdkit && \
    conda install -y nose  \
                     cairo \
                     nomkl \
                     rdkit \
                     uwsgi
RUN pip install boltons jinja2

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
