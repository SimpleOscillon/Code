First, check that your requirements are met by running:
```console
  pip install -r requirements.txt
```
Second, modify the global variables in `lifecycle.py` to the desired values, for example:
```python
############################################################################
# Edit these parameters:
############################################################################

# the values of the frequency to consider:
w_range = np.linspace(0.7, 0.99, 100)
# the Fourier coefficients of the potential. If they do not sum to one,
# another one will be added to satisfy the sum:
coeffs = np.array([0.7, -0.5, 0.5])
# the size of the spatial box:
L = 20.0
# the spatial step size:
dr = 0.01
# number of perturbative harmonics:
N_harmonics = 3
```
and then simply run:
```console
  python lifecycle.py
```


The code will return the corresponding lifetime over this grid of frequencies, in this case:
```console
log10(lifetime)= 6.132809799915554
```
And a plot of the power curve, highlighting the parts where energy is decreasing (blue) or increasing (red):
![image](https://user-images.githubusercontent.com/25751555/112235844-addea080-8bfc-11eb-9a42-6ae580a7ff32.png)

