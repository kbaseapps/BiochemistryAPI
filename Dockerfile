FROM kbase/sdkbase2:python
MAINTAINER KBase Developer2

# -----------------------------------------

RUN apt-get update -y && \
    apt-get install -y build-essential libssl1.1 libssl-dev

RUN conda config --add channels  https://conda.anaconda.org/rdkit && \
    conda install -y uwsgi \
		     nose \
                     cairo \
                     nomkl \
                     rdkit
RUN pip install boltons

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
