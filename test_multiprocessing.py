from multiprocessing import Process, cpu_count, Queue
import time
import numpy as np
def counter(num):
    count=0
    while count<num:
        count+=1
def eat(q,x):
    y=x**2
    q.put(y)
print(cpu_count())
time_0=time.perf_counter()
q=Queue()
a=Process(target=eat,args=(q,np.array([25,30])))
a.start()
b=Process(target=eat,args=(q,np.array([30,50])))
b.start()
#c=Process(target=counter,args=(250000,))
#c.start()
#d=Process(target=counter,args=(250000,))
#d.start()
print('first=',q.get())
print('second=',q.get())
a.join()
b.join()
#c.join()
#d.join()
print('finished in',time.perf_counter()-time_0)
