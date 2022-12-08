import threading
import time
time_0=time.perf_counter()
def eat_breakfast():
    time.sleep(3)
    print('eat_breakfast')

def drink_coffee():
    time.sleep(4)
    print('drink_coffee')

def study():
    time.sleep(5)
    print('study')

x = threading.Thread(target=eat_breakfast,args=None)
x.start()
y = threading.Thread(target=drink_coffee,args=None)
y.start()
z = threading.Thread(target=study,args=None)
z.start()

#eat_breakfast()
#drink_coffee()
#study()

x.join()
y.join()
z.join()

print(threading.active_count())
print(threading.enumerate())
print(time.perf_counter()-time_0)
