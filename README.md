bullard-gellman_integrals
=========================

This project efficiently evaluates Adams Gaunt and Elsasser integrals which are usefull in Geomagnetic calculations. These are of the form

Adams-Gaunt (http://bit.ly/1wiAxzX ) : \int\limits_\Omega Y_{l_3}^{m_3} Y_{l_2}^{m_2} Y_{l_1}^{m_1 *} \sin \theta d\theta d\phi
Elsasser (http://bit.ly/1vU1T0o ): \int\limits_\Omega \left(\frac{\partial Y_{l_2}^{m_2}}{\partial \phi} \frac{\partial Y_{l_3}^{m_3}}{\partial \theta}-\frac{\partial Y_{l_2}^{m_2}}{\partial \theta} \frac{\partial Y_{l_3}^{m_3}}{\partial \phi}\right)Y_{l_1}^{m_1 *}d\theta d\phi 

Where each Y is a spherical harmonic, and each can have a complex conjugate. Warning, use this software at your own risk
