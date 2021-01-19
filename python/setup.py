import configparser
import os.path

import pybind11
import setuptools


def get_metadata(key):
    with open(os.path.join("..", "metadata", key+".txt"), "r", encoding="utf8") as f:
        return f.read().strip()


if __name__ == "__main__":
    metadata = {
        "name": "vect2rast",
        "version": get_metadata("version"),
        "author": get_metadata("author"),
        "author_email": "email",
        "description": get_metadata("description"),
        "url": get_metadata("repository"),
    }

    with open(os.path.join("..", "README.md"), "r") as f:
        metadata["long_description"] = f.read()

    config = configparser.ConfigParser()
    config.read("setup.cfg")
    vect2rast_include_dir = config["vect2rast"].get("include_dir", "")
    vect2rast_library_dir = config["vect2rast"].get("library_dir", "")

    vect2rast = setuptools.Extension(
        "vect2rast.vect2rast",
        include_dirs=[pybind11.get_include(),
                      vect2rast_include_dir],
        sources=[os.path.join("vect2rast",
                              "vect2rast.cpp")],
        libraries=["vect2rast"],
        library_dirs=[vect2rast_library_dir],
    )

    setuptools.setup(
        long_description_content_type="text/markdown",
        packages=setuptools.find_packages(),
        ext_modules=[vect2rast],
        **metadata
    )
