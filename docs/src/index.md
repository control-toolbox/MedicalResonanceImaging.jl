# Optimal control in Medical Resonance Imaging

## Introduction

From the 2010 survey[^1], one of the earliest examples of coherent control of quantum dynamics 
is manipulation of nuclear spin ensembles using radiofrequency (RF) fields[^2]. 
This manipulation is possible due to the Nuclear Magnetic Resonance (NMR) phenomenon[^3] [^4]
which has become a very powerful tool to study matter in a variety of domains from biology 
and chemistry to solid physics and quantum mechanics. Two of the main applications of NMR 
control techniques are maybe high-resolution spectroscopy and Magnetic Resonance Imaging (MRI).

MRI is a medical imaging technique used to produce pictures of the anatomy and investigate 
physiological processes of the body. The general principles involve the interaction of matter 
with electromagnetic fields and are the following. When a sample of matter, in liquid phase, 
is inside a strong and uniform longitudinal magnetic field ($B_0$), the magnetic moments
of the spin-$1/2$ particles align with the direction of this field. When a transverse
radio-frequency magnetic pulse ($B_1$) is then applied, the sample alter its magnetization 
alignment relative to the field and its characteristics (the relaxation times $T_1$ and
$T_2$). These changes in magnetization alignment produce a changing magnetic flux, which 
yields to the emission of an electric signal which is detected and amplified. Then, the image 
is obtained by the analysis of the Fourier transform of this signal. One key point of this 
process is the control of the magnetization vector via the magnetic pulse.
To do that, many different control strategies already exists[^4], but the majority are based 
upon heuristic reasoning[^5] [^6].

![MRI](mri-resources/mri.jpg)

Optimal control algorithms were introduced in NMR to improve control field sequences 
recently[^7] and at the end of the nineties, new methods appeared in optimal control of NMR 
systems both from the analytical and numerical points of view[^8] [^9]. More recently, the 
combination of geometric optimal control based on the Maximum Principle[^10] and related 
numerical algorithms (gradient methods[^11], shooting and continuation methods[^12]) leads 
to sophisticated results about the time-minimal saturation problem (which consists in bringing 
the magnetization vector to zero, *i.e.* the state to the origin) of a single spin[^13]
with applications to the contrast problem in MRI, see[^14] [^15]. They are the basis to 
numerical computations of *robust optimal controls* which take into account inhomogeneities 
contained in the $B_0$ and $B_1$ magnetic fields and which have been validated very recently 
by *in vitro* and *in vivo* experiments, see for instance [^16] [^17].

!!! note "2022 survey"

    This background overview is not up-to-date. We refer to [^18] for more details.

[^1]: C. Brif, R. Chakrabarti & H. Rabitz, *Control of Quantum Phenomena: Past, Present and Future*, New Journal for Physics, **12** (2010), pp.1--68.

[^2]: A. Abragam, *The Principles of Nucelar Magnetism*, Oxford University Press, London (1961).

[^3]: R. R. Ernst, *Principles of Nuclear Magnetic Resonance in one and two dimensions* International Series of Monographs on Chemistry, Oxford University Press, Oxford (1990).

[^4]: M. H. Levitt, *Spin dynamics: Basics of Nuclear Magnetic Resonance*, John Wiley & Sons, New York-London-Sydney (2008).

[^5]: M. A. Berstein, K. F. King & X. J. Zhou, *Handbook of MRI pulse sequences*, Elsevier, Amsterdam (2004).

[^6]: J. N. Rydberg, S. J. Riederer, C. H. Rydberg & C. R. Jack, *Contrast optimization of fluid-attenuated inversion recovery (flair) imaging*, Magn. Reson. Med., **34** (1995), no. 6, 868--877.

[^7]: S. Conolly, D. Nishimura & A. Macovski, *Optimal control solutions to the magnetic resonance selective excitation problem*, IEEE Trans. Med. Imaging, **5** (1986), no. 2, 106--115.

[^8]: N. Khaneja, S. J. Glaser & R. Brockett, *Sub-Riemannian geometry and time optimal control of three spin systems: quantum gates and coherence transfer*, Phys. Rev. A, **65** (2002), no. 3, 032301, pp.11. Errata Phys. Rev. A **68** (2003), 049903; Phys. Rev. A **71** (2005), 039906.

[^9]: T. E. Skinner, T. O. Reiss, B. Luy, N. Khaneja & S. J. Glaser, *Application of optimal control theory to the design of broadband excitation pulses for high-resolution NMR*, J. Magn. Reson., **163** (2003), no. 1, 8--15.

[^10]: L. S. Pontryagin, V. G. Boltyanskii, R. V. Gamkrelidze & E. F. Mishchenko, *The Mathematical Theory of Optimal Processes*, Translated from the Russian by K. N. Trirogoff, edited by L. W. Neustadt, Interscience Publishers John Wiley & Sons, Inc., New York-London, 1962.

[^11]: N. Khaneja, T. Reiss, C. Kehlet, T. Schulte-Herbriiggen & S. J. Glaser, *Optimal control of coupled spin dynamics: design of NMR pulse sequences by gradient ascent algorithms*, J. Magn. Reson., **172** (2005), no. 2, 296--305.

[^12]: O. Cots, *Contrôle optimal géométrique : méthodes homotopiques et applications*, Phd thesis, Institut Mathématiques de Bourgogne, Dijon, France, 2012.

[^13]: M. Lapert, Y. Zhang, M. Braun, S. J. Glaser & D. Sugny, *Singular extremals for the time-optimal control of dissipative spin 1/2 particles*, Phys. Rev. Lett., **104** (2010), 083001.

[^14]: B. Bonnard, M. Chyba & J. Marriott, *Singular Trajectories and the Contrast Imaging Problem in Nuclear Magnetic resonance*, SIAM J. Control Optim., **51** (2013), no. 2, 1325--1349.

[^15]: M. Lapert, Y. Zhang, M. A. Janich, S. J. Glaser & D. Sugny, *Exploring the physical limits of saturation contrast in magnetic resonance imaging*, Sci. Rep., **589** (2012), no. 2.

[^16]: B. Bonnard, O. Cots, S. Glaser, M. Lapert, D. Sugny & Y. Zhang, *Geometric optimal control of the contrast imaging problem in nuclear magnetic resonance, IEEE Trans. Automat. Control, **57** (2012), no. 8, 1957--1969.

[^17]: E. Van Reeth, H. Ratiney, M. Tesch, D. Grenier, O. Beuf, S. J. Glaser &  D. Sugny, *Optimal control design of preparation pulses for contrast optimization in MRI*, J. Magn. Reson., **279** (2017), 39--50.

[^18]: C. P. Koch, U. Boscain, T. Calarco, G. Dirr, S. Filipp, S. J. Glaser, R. Kosloff, S. Montangero, T. Schulte-Herbrüggen, D. Sugny & F. K. Wilhelm, *Quantum optimal control in quantum technologies. strategic report on current status, visions and goals for research in europe*, EPJ Quantum Technology, 9 (2022), p. 60.

## Dependencies

All the numerical simulations to generate this documentation from `MRI.jl` are performed with 
the following packages.

```@example
using Pkg
Pkg.status()
```