# The main function of the fourth edition

What's new since the third version?

1. Add modes whose frequencies are higher than the VIV mode to deal with the frequency multiplier effect





## The modal property 

|      | Order | Frequency (Hz) | Vertical bending                                             |
| ---- | ----- | -------------- | ------------------------------------------------------------ |
| 1    | 2     | 0.093671       | √                                                            |
| 2    | 3     | 0.101788       | √                                                            |
| 3    | 5     | 0.132033       | √                                                            |
| 4    | 6     | 0.177923       | √                                                            |
| 5    | 7     | 0.179536       | √                                                            |
| 6    | 13    | 0.227389       | √                                                            |
| 7    | 20    | 0.271895       | √                                                            |
| 8    | 22    | **0.321445**   | ![image-20240321111548194](Read%20Me.assets/image-20240321111548194.png) |
| 9    | 27    | 0.370955       | √                                                            |
| 10   | 33    | 0.423803       | √                                                            |
| 11   | 39    | 0.476438       | √                                                            |
| 12   | 43    | 0.523472       | √                                                            |
| 13   | 44    | 0.547607       | √                                                            |
| 14   | 46    | 0.587517       | √                                                            |
| 15   | 49    | 0.595777       | √                                                            |
| 16   | 52    | **0.646153**   | ![image-20240321111633673](Read%20Me.assets/image-20240321111633673.png) |
| 17   | 59    | 0.691685       | √                                                            |
| 18   | 61    | 0.738838       | √                                                            |
| 19   | 65    | 0.792543       | √                                                            |
| 20   | 69    | 0.847464       | √                                                            |
| 21   | 71    | 0.900707       | √                                                            |
| 22   | 76    | 0.912536       | √                                                            |
| 23   | 78    | **0.965064**   | ![image-20240321111703082](Read%20Me.assets/image-20240321111703082.png) |
| 24   | 83    | 1.026025       | √                                                            |
| 25   | 87    | 1.089473       | √                                                            |
| 26   | 94    | 1.158499       | √                                                            |
| 27   | 98    | 1.229351       | √                                                            |



## VIV in the frequency domain

### Section 1 

![image-20240321112336064](Read%20Me.assets/image-20240321112336064.png)

### Section 2

![image-20240321112400856](Read%20Me.assets/image-20240321112400856.png)

### Section 3

![image-20240321112434575](Read%20Me.assets/image-20240321112434575.png)

## Apply a bandpass filter [0.63,0.66]

![image-20240321112808539](Read%20Me.assets/image-20240321112808539.png)

#### The front part

![image-20240321113957552](Read%20Me.assets/image-20240321113957552.png)

#### The middle part

![image-20240321112827115](Read%20Me.assets/image-20240321112827115.png)

#### The last part

![image-20240321112938896](Read%20Me.assets/image-20240321112938896.png)

### Apply a bandpass filter [0.95,1]

![image-20240321114100248](Read%20Me.assets/image-20240321114100248.png)

#### The front part

![image-20240321114125655](Read%20Me.assets/image-20240321114125655.png)

#### The middle part

![image-20240321114152228](Read%20Me.assets/image-20240321114152228.png)

### The last part

![image-20240321114247068](Read%20Me.assets/image-20240321114247068.png)

## Apply Kalman filter

![image-20240321144249131](Read%20Me.assets/image-20240321144249131.png)

![image-20240321153752543](Read%20Me.assets/image-20240321153752543.png)
