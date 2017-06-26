# Purpose: generator input file for feap
#  nerbs, elmt01.f
#  elastic

import os
import numpy as np

def  geta(Temperature):
    a = np.zeros(10)
    a[0]=4.124e5*(Temperature-388)
    a[1]=-209.7e6
    a[2]=7.974e8
    a[3]=129.4e7
    a[4]=-1.95e9
    a[5]=-2.5009e9
    a[6]=3.863e10
    a[7]=2.529e10
    a[8]=1.637e10
    a[9]=1.367e10
    return a[np.array([1,2,3,4,5,7,8,9])]


dir_path = os.path.dirname(os.path.realpath(__file__))
Input = open(dir_path+"/Input",'w')

# Set up parameters:
ma=1
Youngs = 1.6e11
poission = 0.32
body_x = 0.0
body_y = 0.0
permittivity = 2000
d31 =-113e-12
d33 = 360e-12
d15 = 162e-12
T=300
coef_a = np.zeros(8)
coef_a = geta(T)
G11 = 1
G12 = 1
G44 = 0.5

# Set up geometry
left_bound = 0.0
right_bound = 1.0
bot_bound = 0.0
top_bound = 5.0

# Set up block number
block_x= 20
block_y= 20


Input.write("feap * * 2D block/ IGA\n")
Input.write("  ndm = 2\n")
Input.write("  ndf = 5\n")
Input.write("\n")
Input.write("INCLude Mesh2d.dat\n") 
Input.write("\n")
Input.write("mate,1\n")
Input.write("user,3\n")
Input.write(str(Youngs)+" "+str(poission)+"\n")
Input.write("1.0\n")
Input.write("1.0\n")
Input.write("3 3\n")  #quadratic
Input.write(str(body_x)+" "+str(body_y)+"\n")
Input.write(str(permittivity*8.85e-12)+"\n")
Input.write(str(d31)+" "+str(d33)+" "+str(d15)+"\n")
Input.write(str(G11)+" "+str(G12)+" "+str(G44)+"\n")
Input.write(str(coef_a).replace('\n', '')[1:-1]+"\n")
Input.write("\n")
Input.write("ebound\n")
Input.write("1 "+str(left_bound) +"  0 0 0 "+"  0 "+" 0"+"\n")
Input.write("1 "+str(right_bound)+"  0 0 0 "+"  0 "+" 0"+"\n")
Input.write("2 "+str(bot_bound)  +"  0 0 0 "+"  0 "+" 1"+"\n")
Input.write("2 "+str(top_bound)  +"  0 0 0 "+"  0 "+" 0"+"\n")
Input.write("\n")
Input.write("cbound"+"\n")
Input.write("node 0 0 0 0 0 1 1"+"\n")
Input.write("node "+str(right_bound)+" 0  0 0 0 0 1"+"\n")
Input.write("\n")
Input.write("end\n")
Input.write("\n")
Input.write("opti\n")
Input.write("\n")
Input.write("batch\n")
Input.write("tang,,1\n")
Input.write("stre,node,1,numnp\n")
Input.write("vtko\n")
Input.write("end\n")

# write mesh2d.dat
mesh = open(dir_path+"/Mesh2d.dat",'w')

#write knots part  0 0 0 ... 1 1 1 

knots_x = block_x + 5
knots_y = block_y + 5
np1 = block_x + 2
np2 = block_y + 2
mesh.write("KNOTS\n")
mesh.write("open 1 " + str(knots_x))
mesh.write(" 0.0 0.0 0.0")
interval=1.0/block_x
word_count=6
for i in range(block_x-1):
    mesh.write(" "+str((i+1)*interval))
    word_count=word_count+1
    if (((word_count) % 16) == 0):
        mesh.write("\n")
for i in range(3):
    mesh.write(" 1.0")
    word_count=word_count+1
    if (((word_count) % 16) == 0):
           mesh.write("\n")
if (((word_count) % 16) != 0):
           mesh.write("\n")
mesh.write("open 2 " + str(knots_y))
mesh.write(" 0.0 0.0 0.0")
interval=1.0/block_y
word_count=6
for i in range(block_y-1):
    mesh.write(" "+str((i+1)*interval))
    word_count=word_count+1
    if (((word_count) % 16) == 0):
        mesh.write("\n")
for i in range(3):
    mesh.write(" 1.0")
    word_count=word_count+1
    if (((word_count) % 16) == 0):
           mesh.write("\n")
mesh.write("\n")

mesh.write("\n")

#write Nurb part  

mesh.write("NURBs\n")
weight=1.0
interval_x=(right_bound-left_bound)/block_x
interval_y=(top_bound-bot_bound)/block_y
count = 1
cp = np.zeros((np1,np2), dtype=np.int)
buff_x = 0.0
for i in range(np1):
    if (i==0 or i==block_x):
        xinc=interval_x/2.0
    else:
        xinc=interval_x
    buff_y = 0.0
    for j in range(np2):
        if (j==0 or j==block_y):
            yinc=interval_y/2.0
        else:
            yinc=interval_y   
        mesh.write(str(count)+" 0 "+str(buff_x)+" "+str(buff_y)+" "+str(weight)+"\n")
        buff_y=buff_y+yinc
        cp[i][j] = count
        count=count+1    
    buff_x=buff_x+xinc
        

#write npatch part
mesh.write("\n\n")
mesh.write("NPATch\n")
mesh.write("  SURFace   "+ str(ma)+" " +str(np1)+" "+str(np2)+" "+" 1 2 " " \n")
#for i in range(np2):
    #for j in range(np1):
      #mesh.write(str(cp[i][j])+" ")
    #mesh.write("\n")


for i in range(np1):
    word_count=0
    mesh.write("     ")
    for j in range(np2):
      mesh.write(str(cp[j][i])+" ")
      word_count=word_count+1
      if (((word_count) % 16) == 0):
           mesh.write("\n")
    if (((word_count) % 16) != 0):
           mesh.write("\n")
    


