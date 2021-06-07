  
FROM debian:unstable
RUN apt-get clean && apt-get update && apt-get install -y cppcheck
COPY . /pipeline/source
RUN cd /pipeline/source && cppcheck --error-exitcode=1 --enable=performance,portability .
