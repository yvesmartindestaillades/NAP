DOCKER_IMAGE := ydmt/dreem_nap
VERSION := $(shell git describe --always --dirty --long)
default:
	echo "See readme"
	pip3 uninstall dreem_nap -y
	pip3 install .
	
init:
	pip install -r requirements.txt

build-image:
	docker build .
		-f ./Dockerfile
		-t $(DOCKER_IMAGE):$(VERSION)

push-image:
	docker push $(DOCKER_IMAGE):$(VERSION)

pin-dependencies:
	pip install -U pip-tools
	pip-compile requirements.in

upgrade-dependencies:
	pip install -U pip pip-tools
	pip-compile -U requirements.in > requirements.txt

push_to_pypi:
	python3 -m build
	twine upload -r pypi dist/*


