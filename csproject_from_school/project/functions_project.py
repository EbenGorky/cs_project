import math
from math import *
#dharmesh
def quantization_of_charge__q(m):
	return("Q=",float(m)*(1.6),"10^-19C\nFormula:Q=±ne")
def number_of_electrons_in_given_charge__n(m):
	return("n=",float(m)/(1.6),"10^19\nFormula:n=Q/e")
def Force_between_two_charges__F(m,n,o):
	return("F=",(float(m)*float(n)/float(o)**2)*9,"10^9FN\nFormula:F=k*q1*q2/r^2")
def distance_between_two_charges__r(m,n,o):
	return("r=",sqrt((float(m)*float(n)*9)*10**9/float(o)),"m\nFormula:r=(q1*q2*k/F)^1/2")
def first_charge__q1(m,n,o):
	return("q1=",float(m)**2*float(o)/(float(n)*9*10**9),"C\nFormula:q1=F*r^2/q2*k")
def second_charge__q2(m,n,o):
	return("q2=",float(n)**2*float(o)/(float(m)*9*10**9),"C\nFormula:q2=F*r^2/q1*k")
def Relation_between_F_and_E__F(m,n):
	return("F=",float(m)*float(n),"N\nFormula:F=qE")
def Electric_field__E(m,n):
	return("E=",float(n)/float(m),"v/m\nFormula:E=F/q")
def charge__q(m,n):
	return("q=",float(n)/float(m),"C\nFormula:q=F/E")
def Electric_field_due_to_a_point_charge__E(m,n):
	return("E=",9*m/n**2,"*10^9v/m\nFormula:E=k*Q/r^2")
def distance_between_charge_and_point__r(m,n):
	return("r=",sqrt((9*10**9)*float(m)/float(n)),"m\nFormula:r=(k*Q/E)^1/2")
def charge__q(m,n):
	return("Q=",(float(n)*float(m)**2)/(9*10**9),"C\nFormula:Q=E*r^2/k")
def Gauss_theorm__phi(m):
	return("ø=",float(m)/8.854,"*10^12\nFormula:ø=q/ε0")
def charge__q(m):
	return("q=",8.854*float(m)*10**-12," C\nFormula:q=ø*ε0")
def Electric_potential_due_to_a_point_charge__V(m,n):
	return("V=",float(m)*9/float(n),"*10^9\nFormula:V=k*q/r")
def charge__q(m,n):
	return("q=",(float(m)*float(n))/(9*10**9),"C\nFormula:q=v*r/k")
def distance__r(m,n):
	return("r=",9*10**9*float(m)/float(n),"m\nFormula:r=k*q/V")
def Electric_potential_due_to_dipole__V(m,n,o):
	return("v=",9*float(m)*cos(n)/(o**2),"*10^9\nFormula:V=k*p*cosθ/r^2")
def distance_between_two_charges__r(m,n,o):
	return("r=",sqrt((9*10**9)*(float(m)*cos(float(n)))/float(o)),"m\nFormula:r=(k*p*cosθ/v)^1/2")
def dipole_moment__p(m,n,o):
	return("p=",float(m)*float(o)**2/(9*10**9*cos(float(n))),"Cm\nFormula:p=v*r^2/k*cosθ")
def angle__theta(m,n,o):
	return("θ=",acos(((float(m)*10**-9)*(float(o)**2))/9*float(n)),"°\nFormula:θ=cos^-1(v*r^2/k*p)")
def Potential_energy_of_a_system_of_two_point_charges__U(m,n,o):
	return("U=",9*float(m)*float(n)/float(o),"*10^9 J\nFormula:U=k*q1*q2/r")
def distance_between_two_point_charges__r(m,n,o):
	return("r=",(9*10**9)*float(m)*float(n)/float(o),"m\nFormula:r=k*q1*q2/U")
def first_charge__q1(m,n,o):
	return("q1=",float(m)*float(o)/(9*10**9*float(n)),"C\nFormula:q1=U*r/k*q2")
def second_charge__q2(m,n,o):
	return("q2=",float(n)*float(o)/(9*10**9*float(m)),"C\nFormula:q2=U*r/k*q1")
def Field_intensity_due_to_infinitely_long_straight_uniformly_charged_wire__E(m,n):
	return("E=",float(m)/(55.603*10**-12*float(n)),"v/m\nFormula:E=λ/2ε0*R")
def Resistance__R(m,n):
	return("R=",float(m)/(55.603*10**-12*float(n)),"m\nR=λ/2πε0*E")
def lamda_value__lamda(m,n):
	return("λ=",float(m)*55.603*10**-12*float(n),"\nFormula:λ=R*2πε0*E")
def Field_intensity_due_to_uniformly_charged_spherical_shell_out__E(m,n):
	return("E=",9*float(m)/(float(n)**2),"*10^9V/m\nFormula:E=k*q/r^2")
def radius_of_shell_including_thickness__r(m,n):
	return("r=",sqrt(((float(m)*9*10**9)))/sqrt(float(n)),"m\nFormula:r=√(q*k/E)")
def charge__q(m,n):
	return("q=",(float(m)*float(n)**2)/9*10**-9,"C\nFormula:q=E*r^2/k")
def Field_intensity_due_to_uniformly_charged_spherical_shell_on__E(m,n):
	return("E=",9*float(m)/(float(n)**2 ),"*10^9V/m\nFormula:E=k*q/R^2")
def radius_of_shell__r(m,n):
	return("R=",sqrt(9*10**9*float(m))/sqrt(float(n)),"m\nFormula:R=√(k*q/E)")
def charge__q(m,n):
	return("q=",float(m)*float(n)**2/(9*10**9),"C\nFormula:q=ER^2/k")
def Field_intensity_due_to_thin_infinite_plane_sheet_of_charge__E(m):
	return("E=",float(m)/8.854,"10^12V/m\nFormula:E=s/ε0")
def conductivity__sigma(m):
	return("s=",float(m)*8.854*10**-12,"ohm^-1\nFormula:E=s/ε0")
def To_find_the_current_in_a_current_carrying_conductor__I(m,n):
	return("I=",float(m)/float(n),"A\nFormula:I=Q/t")
def charge__Q(m,n):
	return("Q=",float(m)*float(n),"C\nFormula:Q=I*t")
def time__t(m,n):
	return("t=",float(m)/float(n),"seconds\nFormula:t=Q/I")
def Ohms_law__V(m,n):
	return("V=",float(m)*float(n),"\nFormula:V=I*R")
def Current__I(m,n):
	return("I=",float(m)/float(n),"A\nFormula:I=V/R")
def Resistance__R(m,n):
	return("R=",float(n)/float(m),"ohm\nFormula:R=V/I")
def Relation_between_R_and_r__R(m,n,o):
	return("R=",float(m)*float(o)/float(n),"ohm\nFormula:R=r*l/A")
def Resistivity__Rho(m,n,o):
	return("r=",float(m)*float(n)/float(o),"ohm metre\nFormula:r=R*A/l")
def Area_of_cross_section__A(m,n,o):
	return("A=",float(m)*float(o)/float(n),"m^2\nFormula:A=rl/A")
def Lenght_of_the_conductor__l(m,n,o):
	return("l=",float(o)*float(n)/float(m),"m\nFormula:l=RA/r")
def Relation_between_R_and_C__C(m):
	return("C=",1/float(m),"ohm^-1\nFormula:C=1/R")
def Resistance__R(m):
	return("R=",1/float(m),"ohm\nFormula:R=1/C")
def Current_density__j(m,n):
	return("j=",float(m)*float(n),"A/m^2\nFormula:j=s*E")
def conductivity__sigma(m,n):
	return("s=",float(m)/float(n),"\nFormula:s=j/E")
def Electric_field__E(m,n):
	return("E=",float(n)/float(m),"V/m\nFormula:E=j/s")
def Electric_power__P(m,n):
	return("P=",float(m)*float(n),"watts\nFormula:P=V*I")
def Voltage__V(m,n):
	return("V=",float(m)/float(n),"\nFormula:V=P/I")
def Current__I(m,n):
	return("I=",float(n)/float(m),"A\nFormula:I=P/I")
def Magnetic_field_due_to_a_straight_conductor_of_infinite_length__B(m,n,o):
	return("B=",12.57*float(m)*float(n)/(float(o)*2*3.14),"*10^-7 T\nFormula:B=µ0*N*I/2*π*r")
def number_of_turns__N(m,n,o):
	return("N=",float(m)*2*3.14*float(o)/(12.57*10**-7*float(n)),"\nFormula:N=B*2*π*r/µ0*I")
def current__I(m,n,o):
	return("I=",float(n)*2*3.14*float(o)/(12.57*10**-7*float(m)),"A\nFormula:I=B*2*π*r/µ0*N")
def radius__r(m,n,o):
	return("r=",(float(n)*12.57*float(m)*10**-7)/(float(o)*2*3.14),"m\nFormula:r=I*µ0*N/B*2*π")
def Force_acting_on_a_charge_particle_in_magnetic_field__F(m,n,o,p):
	return("F=",float(m)*float(n)*float(o)*sin(float(p)),"N\nFormula:F=B*q*v*sinθ")
def Magnetic_field__B(m,n,o,p):
	return("B=",float(m)/(float(n)*float(o)*sin(float(p))),"T\nFormula:B=F/q*v*sinθ")
def charge__q(m,n,o,p):
	return("q=",float(n)/(float(m)*float(o)*sin(float(p))),"C\nFormula:q=F/B*v*sinθ")
def velocity_of_charge_particle__v(m,n,o,p):
	return("v=",float(o)/(float(m)*float(n)*sin(float(p))),"m/s\nFormula:v=F/Bqsinθ")
def angle__theta(m,n,o,p):
	return("θ=",asin(float(p)/float(m)*float(n)*float(o)),"°\nFormula:θ=sin^-1(F/B*q*v)")
def Lorentz_force__F(m,n,o,p):
	return("F=",float(m)*(float(n)+(float(o)*float(p))),"N\nFormula:F=q*[E+(v*B)]")
def charge__q(m,n,o,p):
	return("q=",float(m)/(float(n)+(float(o)*float(p))),"C\nFormula:q=F/(E+(v*B))")
def Electric_field__E(m,n,o,p):
	return("E=",(float(n)/float(m))-(float(o)*float(p)),"V/m\nFormula:E=(F/q)-(v*B)")
def Velocity_of_charge_particle__v(m,n,o,p):
	return("v=",(-float(n)+(float(o)/float(m)))/float(p),"m/s\nFormula:v=(-E+(F/q))/B")
def Magnetic_field__B(m,n,o,p):
	return("B=",(-float(n)+(float(p)/float(m)))/float(o),"T\nFormula:B=(-E+(F/q))/B")
def Force_on_a_current_carrying_conductor_in_magntic_field__F(m,n,o,p):
	return("F=",float(m)*float(n)*float(o)*sin(float(p)),"N\nFormula:F=B*I*L*sinθ")
def Magnetic_field__B(m,n,o,p):
	return("B=",float(m)/(float(n)*float(o)*sin(float(p))),"T\nFormula:B=F/(I*L*sinθ)")
def current__I(m,n,o,p):
	return("I=",float(n)/(float(m)*float(o)*sin(float(p))),"A\nFormula:I=F/(B*L*sinθ)")
def Lenght__L(m,n,o,p):
	return("L=",float(o)/(float(m)*float(n)*sin(p)),"m\nFormula:L=F/(B*I*sinθ)")
def angle__theta(m,n,o,p):
	return("θ=",asin(float(p)/(float(m)*float(n)*float(o))),"°\nFormula:θ=sin^-1(F/(B*I*L))")
def Torque_on_a_straight_conductor_in_magnetic_field__tau(m,n,o,p,q):
	return("τ=",float(m)*float(n)*float(o)*float(p)*sin(float(q)),"Nm\nFormula:τ=B*I*N*A*sinθ")
def magnetic_field__B(m,n,o,p,q):
	return("B=",float(m)/(float(n)*float(o)*float(p)*sin(float(q))),"T\nFormula:B=τ/(I*N*A*sinθ)")
def current__I(m,n,o,p,q):
	return("I=",float(n)/(float(m)*float(o)*float(p)*sin(float(q))),"A\nFormula:I=τ/(B*N*A*sinθ)")
def number_of_turns__N(m,n,o,p,q):
	return("N=",float(o)/(float(m)*float(n)*float(p)*sin(float(q))),"\nFormula:N=τ/(B*I*A*sinθ")
def Area_A(m,n,o,p,q):
	return("A=",float(p)/(float(m)*float(n)*float(o)*sin(float(q))),"m^2\nFormula:A=τ/(B*I*N*sinθ)")
def angle__theta(m,n,o,p,q):
	return("θ=",asin(float(q)/(float(m)*float(n)*float(o)*float(p))),"°\nFormula:θ=sin^-1(τ/(B*I*N*A))")
def coversion_of_galvanometer_into_voltmeter(m,n,o):
	return("R=",(float(m)/float(n))-float(o),"ohm\nFormula:R=(V/ig)-G")
def voltage__V(m,n,o):
	return("V=",float(n)*(float(m)+float(o)),"\nFormula:")

#eben
def Maximum_amplitude_of_wave_interference__a(a1,a2,fi):
    A = (a1**2 + a2**2 + 2*a1*a2*cos(fi))**(1/2)
    return ('A=',A,':a1**2 + a2**2 + 2*a1*a2*cos(fi)')
def focus_with_given_centure_of_curvature__f(r):
	return ('f=',r/2,':f=r/2')
def focul_length_of_given_distance_between_object_and_image__f(u,v):
	return ('f=',u*v/(u+v),':f=u*v/(u+v)')
def magnifiaction_of_image_in_given_hight_of_object_and_image__m(h,H):
	return ('m=',h/H,':m=h/H')
def magnification_of_image_in_given_distance_between_object_and_image__m(u,v):
	return ('m=',u/v,':m=u/v')
def refractive_index_of_new_medium_with_respect_to_old_medium_snell_law__n(i,r):
	return ('n=',sin(i)/sin(r),':n=sin(i)/sin(r)')
def find_refractive_index_of_medium_one_to_tow_with_tow_to_one__n(n):
	return ('n=',1/n,':n=1/n')
def find_velocity_of_light_in_second_medium_with_i_r_and_v_in_first_medium__v(v,i,r):
	return ('v =', v*sin(r)/sin(i),':v = v*sin(r)/sin(i)')
def critical_angle_in_TIR__i(a,s,i,n):
	return ('i =', asin(n),':i = asin(n)')
def radius_of_curvature_in_given_refractive_index_of_two_medium__r(N,n,u,v):
	return ('r=',(N-n)*u*v/(u*N-v*n),':r=(N-n)*u*v/(u*N-v*n)')
def focal_length_of_lens_in_given_distance_of_oject_and_image__f(u,v):
	return ('f=',u*v/(u-v),':f=u*v/(u-v)')
def distance_of_image_through_lens_in_given_focal_length_and_distance_of_object__v(u,f):
	return ('v=',u*f/(u+f),':v=u*f/(u+f)')
def distance_of_object_through_lens_in_given_focal_length_and_distance_of_image__u(v,f):
	return ('u=',v*f/(f-v),':u=v*f/(f-v)')
def power_of_lens_in_given_focal_length__p(f):
	return ('p=',1/f,':p=1/f')
def power_of_mirror_in_given_focal_length__p(f):
	return ('p=',(-1)/f,':p=(-1)/f')
def focal_length_of_lens_in_given_power__f(p):
	return ('f=',1/p,':f=1/p')
def focal_length_of_mirror_in_given_power__f(p):
	return ('f=',1/p,':f=1/p')
def combinaction_of_two_lens__f(f,F):
	return ('f=',f*F/(f+F),':f=f*F/(f+F)')
def equivalence_magnification_of_two_lens__m(m,M):
	return ('m=',m*M,':m=m*M')
def total_deviation_of_light_in_through_prism__d(i,e,a):
	return ('d=',i+e-a,':d=i+e-a')
def refractive_index_of_prism_in_condisiton_of_light_travelling_from_medium_one_to_two__n(a,d):
	return ('n=',sin((a+d)/2)/sin(a/2),':n=sin((a+d)/2)/sin(a/2)')
def magnification_of_microscope__m(d,f):
	return ('m=',1+d/f,':m=1+d/f')
def magnification_of_compound_microscope__m(l,f):
	return ('m=',l/f,':m=l/f')
def angular_magnification_of_compound_microscope__m(d,f):
	return ('m=',1+d/f,':m=1+d/f')
def wave_length_of_light_in_medium_two_travelling_form_medium_one_to_two__L(V,v,l):
	return ('L=',V/v*l,':L=V/v*l')
def apparent_frequency_dopular_effect_in_given_velocites_of_observer_V_and_source_v__f(c,V,v,f):
	return ('f=',(c+V)/(c-v)*f,':f=(c+V)/(c-v)*f')
def distance_of_fringe_of_a_given_light__x(n,l,D,d):
	return ('x=',n*l*D/d,':x=n*l*D/d')
def resolving_power_of_lens__a(l,t):
	return ('a=',0.61*l/t,':a=0.61*l/t')
def wave_length_of_polarisation_of_light__l(k):
	return ('l=',2*pi/k,':l=2*pi/k')
def intensity_of_polarised_light_in_given_maximum_intensity_of_I__i(I,t):
	return ('i=',I(cos(t))**2,':i=I(cos(t))**2')
def kinetic_energy_of_electron_in_a_incident_of_light__k(v):
	return ('k=',1.9*10**(-19)*v,':k=1.9*10**(-19)*v')
def wave_length_of_partcle__l(m,v):
	return ('l=',6.626*10**(-34)/(m*v),':l=6.626*10**(-34)/(m*v)')
def deBroglie_wave_lenght__l(v):
	return ('l=',1.227/v**0.5,':l=1.227/v**0.5')
def total_energy_of_electron__e(r):
	return ('e=',((-1)*1.9*10**(-19))**2/(8*pi*(8.854187817*10**(-12))*r),':e=((-1)*1.9*10**(-19))**2/(8*pi*(8.854187817*10**(-12))*r)')

def frequency_of_light_emighted_when_electron_travel_form_starting_N_to_ending_n_orbit__f(n,N):
	return ('f=',(1.097*10**7)*(1/n**2-1/N**2),':f=(1.097*10**7)*(1/n**2-1/N**2)')
def velocity_of_electron_in_nth_orbit__v(n):
	return ('v=',(1.9*10**(-19))**2/(n*2*8.854187817*10**(-12)*1.097*10**7),':v=(1.9*10**(-19))**2/(n*2*8.854187817*10**(-12)*1.097*10**7)')
def radius_of_orbit_at_nth_orbit__r(n):
	return ('r=',n**2*(1.097*10**7)**2*8.854187817*10**(-12)/(pi*(9.1093837*10**(-31))*(1.9*10**(-19))**2),':r=n**2*(1.097*10**7)**2*8.854187817*10**(-12)/(pi*(9.1093837*10**(-31))*(1.9*10**(-19))**2)')
def total_energy_of_electron_orbiting_at_nth_orbit__e(n):
	return ('e=',-2.18*10*(-19)/n**2,':e=-2.18*10*(-19)/n**2')
def energy_of_nuclear_binding_energy__e(m):
	return ('e=',m*(3*10**8)**2,':e=m*(3*10**8)**2')
def radioactivity_active_decay_at_given_time_t__n(N,l,t):
	return ('n=',N*e**(-l*t),':n=N*e**(-l*t)')
