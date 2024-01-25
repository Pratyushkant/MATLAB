
After downloading the MATLAB zip file form the official website, open the directory where the zip file is located. Unzip it using the following command:

~~~bash
mint@Dell-G:~/Downloads$ unzip matlab_R2023b_glnxa64.zip 
~~~


After that type in the command

~~~bash
sudo ./install
~~~

Enter your root password and the MATLAB installer will open. Follow the on-screen instructions to download and install the MATLAB.

To activate the MATLAB, use the following commands:

~~~bash
mint@Dell-G:~$ cd /usr/local/MATLAB/R2023b/bin/glnxa64
mint@Dell-G:/usr/local/MATLAB/R2023b/bin/glnxa64$ sudo ./MathWorksProductAuthorizer
[sudo] password for mavel_math:                                 
mint@Dell-G:/usr/local/MATLAB/R2023b/bin/glnxa64$ 
~~~

To access the MATLAB, type the command

~~~bash
mint@Dell-G:/usr/local/MATLAB/R2023b/bin$ ./matlab
~~~

in the terminal in the MATLAB/bin directory. Usually, you have to access the MATLAB only via the terminal, there won't be any graphical interface for doing so.
