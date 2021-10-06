# Python code to draw the required plot

from matplotlib.markers import MarkerStyle
import pandas as pd
from scipy.optimize import curve_fit
from scipy.odr import Model, RealData, ODR
import matplotlib.pyplot as plt
import numpy as np

def powerlaw(x,a,c):
    f=a*x**c
    return f

def angular_to_solar(RelativeRadius):
    return Distance*StellarRadius*RelativeRadius*0.215/2

def solar_to_angular(RealSize):
    return 2*RealSize/(0.215*StellarRadius*Distance)

def f(B, x):
    return B[0]*x**B[1]

# Star="RDor"
# Star="WHya"
Star="Both"
global Distance
global StellarRadius

if Star=="RDor" or Star=='Both':
    Distance=55
    StellarRadius=55

    datafile="Star_Tables.csv/R_Dor.csv"

    df = pd.read_csv(datafile)

    """ ODR fitting: """
    model = Model(f)
    mydata = RealData(x=df["r/R"], y=df["Brightness Temperature [K]"], sx=df["Relative Size Uncertainty"], sy= df["Bright.Temp. Uncertainty [K]"])
    myodr = ODR(mydata, model, beta0=[3000., -2.])
    myoutput = myodr.run()
    pars1=myoutput.beta
    cov1=myoutput.cov_beta
    myoutput.pprint()
    print("\n ODR Betas:", myoutput.beta)
    '''----------------'''

    x_space=np.linspace(1,1.5, num=300)
    pars, cov = curve_fit(powerlaw, df["r/R"], df["Brightness Temperature [K]"], sigma=df["Bright.Temp. Uncertainty [K]"],p0=[2000,-2])

    fig1 = plt.Figure(figsize=(6, 6))
    ax1 = fig1.subplots()
    
    ax1.plot(x_space, powerlaw(x_space, *pars1),    c='r', zorder=1, label="Orthogonal Distance Regression (ODR)", alpha=0.6, linewidth=0.9)
    ax1.plot(x_space, powerlaw(x_space, *pars),     c='g', zorder=1, label="Ordinary Least Squares (OLS)", alpha=0.6, linewidth=0.9)
    ax1.errorbar(df["r/R"], df["Brightness Temperature [K]"], xerr=df["Relative Size Uncertainty"] ,yerr=df["Bright.Temp. Uncertainty [K]"], fmt='o', ecolor='k', ls = "None", color="None", elinewidth=0.5, capsize=2, zorder=2)
    ax1.scatter(df["r/R"], df["Brightness Temperature [K]"], c='k', s=6, zorder=3)
    
    textstr = '\n'.join((
        "\nOLS fit power-law:",
        r'$T_b= (%.0f \pm %.0f)_K \times {\left(\frac{r}{R_*}\right)}^{(%.2f \pm %.2f)}$' % (pars[0], np.sqrt(cov[0][0]), pars[1], np.sqrt(cov[1][1])),
        r"ODR fit power-law:",
        r'$T_b= (%.0f \pm %.0f)_K \times {\left(\frac{r}{R_*}\right)}^{(%.2f \pm %.2f)}$' % (pars1[0], np.sqrt(cov1[0][0]), pars1[1], np.sqrt(cov1[1][1])),
    ))
    props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor="gray")
    ax1.text(0.48, 0.69, textstr, transform=ax1.transAxes, fontsize=10, horizontalalignment='left',\
            verticalalignment='bottom', bbox=props)
    ax1.text(0.02, 0.98, "R_Dor", fontsize='xx-large', transform=ax1.transAxes, horizontalalignment='left',\
            verticalalignment='top')

    ax1.legend(fontsize=10)

    ax1.tick_params(axis="y",direction="in")
    ax1.tick_params(axis="x",direction="in")

    ax1.set_xlabel(r"Relative radius $\left( r/R_*\right)$",fontsize=12)
    ax1.set_ylabel(r"Brightness Temperature [K]",fontsize=12)
    ax1.set_xlim((0.99,1.31))
    ax1.set_ylim((1550,2650))

    secax= ax1.secondary_xaxis('top', functions=(angular_to_solar, solar_to_angular))
    secax.tick_params(direction="in")
    secax.set_xlabel(r"$R/R_\odot$", fontsize=12)

    fig1.tight_layout()

    fig1.savefig("RDor_TempSize.pdf")
    plt.close()

if Star=="WHya" or Star=="Both":
    Distance=98
    StellarRadius=40

    datafile="Star_Tables.csv/W_Hya.csv"

    df = pd.read_csv(datafile)

    """ ODR fitting: """
    model = Model(f)
    mydata = RealData(x=df["r/R"], y=df["Brightness Temperature [K]"], sx=df["Relative Size Uncertainty"], sy= df["Bright.Temp. Uncertainty [K]"])
    myodr = ODR(mydata, model, beta0=[3000., -2.])
    myoutput = myodr.run()
    pars1=myoutput.beta
    cov1=myoutput.cov_beta
    myoutput.pprint()
    print("\n ODR Betas:", myoutput.beta)
    '''----------------'''

    x_space=np.linspace(0.9,1.5, num=300)
    pars, cov = curve_fit(powerlaw, df["r/R"], df["Brightness Temperature [K]"], sigma=df["Bright.Temp. Uncertainty [K]"],p0=[2000,-2])

    fig1 = plt.Figure(figsize=(6, 6))
    ax1 = fig1.subplots()
    
    ax1.plot(x_space, powerlaw(x_space, *pars1),    c='r', zorder=1, label="Orthogonal Distance Regression (ODR)", alpha=0.6, linewidth=0.9)
    ax1.plot(x_space, powerlaw(x_space, *pars),     c='g', zorder=1, label="Ordinary Least Squares (OLS)", alpha=0.6, linewidth=0.9)
    ax1.errorbar(df["r/R"], df["Brightness Temperature [K]"], xerr=df["Relative Size Uncertainty"] ,yerr=df["Bright.Temp. Uncertainty [K]"], fmt='o', ecolor='k', ls = "None", color="None", elinewidth=0.5, capsize=2, zorder=2)
    ax1.scatter(df["r/R"], df["Brightness Temperature [K]"], c='k', s=6, zorder=3)
    
    textstr = '\n'.join((
        "\nOLS fit power-law:",
        r'$T_b= (%.0f \pm %.0f)_K \times {\left(\frac{r}{R_*}\right)}^{(%.2f \pm %.2f)}$' % (pars[0], np.sqrt(cov[0][0]), pars[1], np.sqrt(cov[1][1])),
        r"ODR fit power-law:",
        r'$T_b= (%.0f \pm %.0f)_K \times {\left(\frac{r}{R_*}\right)}^{(%.2f \pm %.2f)}$' % (pars1[0], np.sqrt(cov1[0][0]), pars1[1], np.sqrt(cov1[1][1])),
    ))
    props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor="gray")
    ax1.text(0.482, 0.69, textstr, transform=ax1.transAxes, fontsize=10, horizontalalignment='left',\
            verticalalignment='bottom', bbox=props)
    ax1.text(0.02, 0.98, "W_Hya", fontsize='xx-large', transform=ax1.transAxes, horizontalalignment='left',\
            verticalalignment='top')

    ax1.legend(fontsize=10)

    ax1.set_xlabel(r"Relative radius $\left( r/R_*\right)$",fontsize=12)
    ax1.set_ylabel(r"Brightness Temperature [K]",fontsize=12)
    ax1.set_xlim((1.092,1.41))
    ax1.set_ylim((2150,3250))

    ax1.tick_params(axis="x",direction="in")
    ax1.tick_params(axis="y",direction="in")

    secax= ax1.secondary_xaxis('top', functions=(angular_to_solar, solar_to_angular))
    secax.tick_params(direction="in")
    secax.set_xlabel(r"$R/R_\odot$", fontsize=12)

    fig1.tight_layout()

    fig1.savefig("WHya_TempSize.pdf")
    plt.close()

