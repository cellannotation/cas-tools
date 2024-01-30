# How to run notebooks

Create a virtual environment:

```
cd notebooks
python3 -m venv venv
source venv/bin/activate

pip install ipykernel
python3 -m ipykernel install --user --name=venv
```

Install and run Jupyter notebook:
```
pip install notebook
jupyter notebook
```