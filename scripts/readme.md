## LEBC Results

This readme should have the results of the code in this github repo for the following particles
- [Single_Sphere](#singlesphere)
- [Rod2](#rod2)
- [Rod4](#rod4)
- [Rod6](#rod6)
- [Multisphere13](#multisphere13)
- [Mixture of Multi and Single Sphere](#mixture)

(todo: add images of these particles, it'll look nice)


<a name="singlesphere"></a>
# Single Sphere

### Stress vs Volume Fraction

![]() <img src="single_sphere/xy_single_sphere.png"  width="60%">

![]() <img src="single_sphere/yy_single_sphere.png"  width="60%">

### Stress vs Time Results

| Volume Fraction | Stress vs Time |
| ------------- | ------------- |
| 0.05 | []() <img src="single_sphere/stress_single_sphere_0.png"  width="50%">  |
| 0.13 | []() <img src="single_sphere/stress_single_sphere_1.png"  width="50%">  |
| 0.23 | []() <img src="single_sphere/stress_single_sphere_2.png"  width="50%">  |
| 0.32 | []() <img src="single_sphere/stress_single_sphere_3.png"  width="50%">  |
| 0.41 | []() <img src="single_sphere/stress_single_sphere_4.png"  width="50%">  |
| 0.50 | []() <img src="single_sphere/stress_single_sphere_5.png"  width="50%">  |


<a name="rod2"></a>
# Rod2

### Stress vs Volume Fraction

![]() <img src="rod2/xy_rod2.png"  width="60%">

![]() <img src="rod2/yy_rod2.png"  width="60%">

### Stress vs Time Results


| Volume Fraction | Stress vs Time |
| ------------- | ------------- |
| 0.025| []() <img src="rod2/stress_rod2_0.png"  width="50%">  |
| 0.05 | []() <img src="rod2/stress_rod2_1.png"  width="50%">  |
| 0.1  | []() <img src="rod2/stress_rod2_2.png"  width="50%">  |
| 0.2  | []() <img src="rod2/stress_rod2_3.png"  width="50%">  |
| 0.3  | []() <img src="rod2/stress_rod2_4.png"  width="50%">  |
| 0.4  | []() <img src="rod2/stress_rod2_5.png"  width="50%">  |
| 0.45 | []() <img src="rod2/stress_rod2_6.png"  width="50%">  |
| 0.50 | []() <img src="rod2/stress_rod2_7.png"  width="50%">  |

<a name="rod4"></a>
# rod4

### Stress vs Volume Fraction

![]() <img src="rod4/xy_rod4.png"  width="60%">

![]() <img src="rod4/yy_rod4.png"  width="60%">

### Stress vs Time Results


| Volume Fraction | Stress vs Time |
| ------------- | ------------- |
| 0.025| []() <img src="rod4/stress_rod4_0.png"  width="50%">  |
| 0.05 | []() <img src="rod4/stress_rod4_1.png"  width="50%">  |
| 0.1  | []() <img src="rod4/stress_rod4_2.png"  width="50%">  |
| 0.2  | []() <img src="rod4/stress_rod4_3.png"  width="50%">  |
| 0.3  | []() <img src="rod4/stress_rod4_4.png"  width="50%">  |
| 0.4  | []() <img src="rod4/stress_rod4_5.png"  width="50%">  |
| 0.45 | []() <img src="rod4/stress_rod4_6.png"  width="50%">  |
| 0.50 | []() <img src="rod4/stress_rod4_7.png"  width="50%">  |



<a name="rod2"></a>
# rod6

### Stress vs Volume Fraction

![]() <img src="rod6/xy_rod6.png"  width="60%">

![]() <img src="rod6/yy_rod6.png"  width="60%">

### Stress vs Time Results


| Volume Fraction | Stress vs Time |
| ------------- | ------------- |
| 0.025| []() <img src="rod6/stress_rod6_0.png"  width="50%">  |
| 0.05 | []() <img src="rod6/stress_rod6_1.png"  width="50%">  |
| 0.1  | []() <img src="rod6/stress_rod6_2.png"  width="50%">  |
| 0.2  | []() <img src="rod6/stress_rod6_3.png"  width="50%">  |
| 0.3  | []() <img src="rod6/stress_rod6_4.png"  width="50%">  |
| 0.4  | []() <img src="rod6/stress_rod6_5.png"  width="50%">  |
| 0.45 | []() <img src="rod6/stress_rod6_6.png"  width="50%">  |
| 0.50 | []() <img src="rod6/stress_rod6_7.png"  width="50%">  |

<a name="multisphere13"></a>
# Multisphere13
ran with MPI 4 4 2

### Stress vs Volume Fraction

![]() <img src="multisphere_13/xy_multisphere_13.png"  width="60%">

![]() <img src="multisphere_13/yy_multisphere_13.png"  width="60%">

### Stress vs Time Results

| Volume Fraction | Stress vs Time |
| ------------- | ------------- |
| 0.05 | []() <img src="multisphere_13/stress_multisphere_13_0.png"  width="50%">  |
| 0.13 | []() <img src="multisphere_13/stress_multisphere_13_1.png"  width="50%">  |
| 0.23 | []() <img src="multisphere_13/stress_multisphere_13_2.png"  width="50%">  |
| 0.32 | []() <img src="multisphere_13/stress_multisphere_13_3.png"  width="50%">  |
| 0.41 | []() <img src="multisphere_13/stress_multisphere_13_4.png"  width="50%">  |
| 0.50 | []() <img src="multisphere_13/stress_multisphere_13_5.png"  width="50%">  |




<a name="mixture"></a>
# Multisphere13 with Same Size Single Sphere
ran with MPI 4 4 2
> Warning: bug in lebc.py for mixtures is causing particles to not get inserted at high volume fractions. This test needs reran after lebc.py is fixed for mixtures.

### Stress vs Volume Fraction

![]() <img src="multisphere_13_and_single_sphere/xy_multi_and_single.png"  width="60%">

![]() <img src="multisphere_13_and_single_sphere/yy_multi_and_single.png"  width="60%">

### Stress vs Time Results

| Volume Fraction | Stress vs Time |
| ------------- | ------------- |
| 0.05 | []() <img src="multisphere_13_and_single_sphere/stress_multi_and_single_0.png"  width="50%">  |
| 0.13 | []() <img src="multisphere_13_and_single_sphere/stress_multi_and_single_1.png"  width="50%">  |
| 0.23 | []() <img src="multisphere_13_and_single_sphere/stress_multi_and_single_2.png"  width="50%">  |
| 0.32 | []() <img src="multisphere_13_and_single_sphere/stress_multi_and_single_3.png"  width="50%">  |
| 0.41 | []() <img src="multisphere_13_and_single_sphere/stress_multi_and_single_4.png"  width="50%">  |
| 0.50 | []() <img src="multisphere_13_and_single_sphere/stress_multi_and_single_5.png"  width="50%">  |

