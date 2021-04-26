#Requires Python >= 3.7 
import time
import subprocess
import hashlib
def cpu_info(): 
  return subprocess.check_output("lscpu", shell=True).strip().decode()
def git_hash(): 
  return subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()

start = time.perf_counter_ns()
time.sleep(1)
end = time.perf_counter_ns()
print("Hello World")
print("Physical cores:", psutil.cpu_count(logical=False))
print(end - start)
print(git_hash())
print(hashlib.sha1(cpu_info().encode()).hexdigest())
print(platform.processor())
print(cpu_info())
#list all files in dir
#Run and time each thing in the directory
#print into csv 
