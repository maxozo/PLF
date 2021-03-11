import sched, time
import os
s = sched.scheduler(time.time, time.sleep)
count=0
def do_something(sc): 
    print(f"Doing stuff...")
    global count
    count+=1
    f= open(f"guru{count}.txt","w+")
    # do your stuff
    s.enter(20, 1, do_something, (sc,))

s.enter(20, 1, do_something, (s,))
s.run()