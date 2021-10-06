# Python code to draw the required plot

from matplotlib.markers import MarkerStyle
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import numpy as np

def powerlaw(x,a,c):
    f=a*x**c
    return f

def angular_to_solar(RelativeRadius):
    return Distance*StellarRadius*RelativeRadius*0.215/2

def solar_to_angular(RealSize):
    return 2*RealSize/(0.215*StellarRadius*Distance)


def freq_to_wavelength(freq):
    return 299.8/freq

def wavelength_to_freq(wavelength):
    return 299.8/wavelength

# Star="RDor"
# Star="WHya"
Star="Both"
global Distance
global StellarRadius

if Star=="WHya" or Star=='Both':
    Distance=98
    StellarRadius=40

    datafile="Star_Tables.csv/W_Hya.csv"

    df = pd.read_csv(datafile)

    x_space=np.linspace(50,450, num=150)
    pars, cov = curve_fit(powerlaw, df["Frequency [GHz]"], df["r/R"], sigma=df["Relative Size Uncertainty"], p0=[0,0])
    print("cov=", cov)

    fig1 = plt.Figure(figsize=(6, 6))
    ax1 = fig1.subplots()

    ax1.plot(x_space, powerlaw(x_space, *pars), c='g', zorder=1, label="best fit", alpha=0.8)
    ax1.errorbar(df["Frequency [GHz]"], df["r/R"], xerr=2 ,yerr=df["Relative Size Uncertainty"], fmt='o', ecolor='k', ls = "None", color="None", elinewidth=0.5, capsize=2, zorder=2)
    ax1.scatter(df["Frequency [GHz]"], df["r/R"], c='k', s=6, zorder=3)
    textstr = '\n'.join((
        "Best fitted power-law:\n"
        r'$r/R_*=(%.3f \pm %.3f) \times \nu^{(%.3f \pm %.3f)}$' % (pars[0], np.sqrt(cov[0][0]), pars[1], np.sqrt(cov[1][1])),
    ))
    props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='gray')
    ax1.text(0.295, 0.85, textstr, transform=ax1.transAxes, fontsize=12, horizontalalignment='left',\
            verticalalignment='center', bbox=props)
    ax1.text(0.06, 0.97, "W_Hya", fontsize='xx-large', transform=ax1.transAxes, horizontalalignment='left',\
            verticalalignment='top')
    
    ax1.legend(fontsize=12)
    ax1.set_xlabel(r"Frequency [GHz] $(\nu)$", fontsize=12)
    ax1.set_ylabel(r"Relative radius $\left( r/R_*\right)$", fontsize=12)
    ax1.set_xlim((90,430))
    ax1.set_ylim((1.17, 1.37))

    ax1.tick_params(axis="x",direction="in")
    ax1.tick_params(axis="y",direction="in")

    print(angular_to_solar(df["r/R"]))
    
    secax= ax1.secondary_yaxis('right', functions=(angular_to_solar, solar_to_angular))
    secax.tick_params(direction="in")
    secax.set_ylabel(r"$~R/R_\odot$", fontsize=12, rotation=270)

    secax2= ax1.secondary_xaxis('top', functions=(freq_to_wavelength, wavelength_to_freq))
    secax2.tick_params(direction="in")
    secax2.set_xlabel(r"$\lambda~\left[ mm \right]$", fontsize=12)


    # ax1.set_title(r"$\lambda~\left[ mm \right]$", fontsize=10)
    fig1.savefig("WHya_SizeFreq.pdf")
    plt.close()

if Star=="RDor" or Star=='Both':
    Distance=55
    StellarRadius=55

    datafile="Star_Tables.csv/R_Dor.csv"

    df = pd.read_csv(datafile)

    x_space=np.linspace(50,500, num=150)
    pars, cov = curve_fit(powerlaw, df["Frequency [GHz]"], df["r/R"], sigma=df["Relative Size Uncertainty"], p0=[0,0])
    print("cov=", cov)

    fig1 = plt.Figure(figsize=(6, 6))
    ax1 = fig1.subplots()

    ax1.plot(x_space, powerlaw(x_space, *pars), c='g', zorder=1, label="best fit", alpha=0.8)
    ax1.errorbar(df["Frequency [GHz]"], df["r/R"], xerr=2 ,yerr=df["Relative Size Uncertainty"], fmt='o', ecolor='k', ls = "None", color="None", elinewidth=0.5, capsize=2, zorder=2)
    ax1.scatter(df["Frequency [GHz]"], df["r/R"], c='k', s=6, zorder=3)
    textstr = '\n'.join((
        "Best fitted power-law:\n"
        r'$r/R_*=(%.3f \pm %.3f) \times \nu^{(%.3f \pm %.3f)}$' % (pars[0], np.sqrt(cov[0][0]), pars[1], np.sqrt(cov[1][1])),
    ))
    props = dict(boxstyle='round', facecolor='white', alpha=0.5, edgecolor='gray')
    ax1.text(0.295, 0.85, textstr, transform=ax1.transAxes, fontsize=12, horizontalalignment='left',\
            verticalalignment='center', bbox=props)
    ax1.text(0.03, 0.97, "R_Dor", fontsize='xx-large', transform=ax1.transAxes, horizontalalignment='left',\
            verticalalignment='top')
    
    ax1.legend(fontsize=12)
    ax1.set_xlabel(r"Frequency [GHz] $(\nu)$", fontsize=12)
    ax1.set_ylabel(r"Relative radius $\left( r/R_*\right)$", fontsize=12)
    ax1.set_xlim((90,430))
    ax1.set_ylim((1.04, 1.249))

    ax1.tick_params(axis="x",direction="in")
    ax1.tick_params(axis="y",direction="in")

    print(angular_to_solar(df["r/R"]))
    
    secax= ax1.secondary_yaxis('right', functions=(angular_to_solar, solar_to_angular))
    secax.tick_params(direction="in")
    secax.set_ylabel(r"$~~R/R_\odot$", fontsize=12, rotation=270)

    secax2= ax1.secondary_xaxis('top', functions=(freq_to_wavelength, wavelength_to_freq))
    secax2.tick_params(direction="in")
    secax2.set_xlabel(r"$\lambda~\left[ mm \right]$", fontsize=12)

    # ax1.set_title(r"$\lambda~\left[ mm \right]$", fontsize=10)
    fig1.savefig("RDor_SizeFreq.pdf")
    plt.close()