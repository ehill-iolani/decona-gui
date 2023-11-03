# decona-gui
A shiny app for running the decona pipeline through a GUI

## Installation
Pull the github repository
```
git clone 
```
Build the docker image
```
docker build -t decona-gui:latest .
```
Run the docker image as a container
```
docker run --name=decona-gui --rm -d -p 8080:3838 decona-gui:latest
```
