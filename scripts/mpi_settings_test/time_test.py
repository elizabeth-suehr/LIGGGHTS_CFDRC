import numpy as np
import matplotlib.pyplot as plt

def get_data(ms,mpi,vf_id):
    ms_file_base = "log.shear.multisphere"
    s_file_base = "log.shear.s_sphere_"
    end_base = ".liggghts"
    names = []
    loop_times = []
    particles = []
    for i, mpi_str in enumerate(mpi):
        count = -1
        loop_times.append([])
        particles.append([])
        for v in vf_id:
            # single spheres
            # print("single sphere: " + str(i) + " " + str(count))
            count +=1
            loop_times[i].append([])
            particles[i].append([])
            filename = s_file_base + mpi_str + str(v) + end_base
            names.append("s_" + str(v))
            file = open(filename,"r")

            if file:
                for line in file.readlines():
                    if "Loop time of" in line:
                        # print(line)
                        value = float(line.split()[3]) / float(line.split()[8])
                        loop_times[i][count].append(value)
                    if "INFO: Particle insertion ins1: inserted" in line:
                        # print(line)
                        value = float(line.split()[5])
                        particles[i][count].append(value)
            else:
                count -=1
                loop_times[i].pop()
                particles[i].pop()
            
            for m in ms:
                # multi spheres
                # print("multi sphere: " + str(i) + " " + str(count))

                count +=1
                loop_times[i].append([])
                particles[i].append([])
                filename = ms_file_base + m + mpi_str + str(v) + end_base
                names.append("m" + m + str(v))
                file = open(filename,"r")
                if file:
                    for line in file.readlines():
                        if "Loop time of" in line:
                            # print(line)
                            value = float(line.split()[3]) / float(line.split()[8])
                            loop_times[i][count].append(value)
                        if "INFO: Particle insertion ins1: inserted" in line:
                            # print(line)
                            value = float(line.split()[5])
                            particles[i][count].append(value)
                else:
                    count -=1
                    loop_times[i].pop()
                    particles[i].pop()
    return names, loop_times, particles







ms = ["13_","43_","93_","163_"]
mpi = ["MPI221_","MPI321_","MPI332_","MPI442_"]
vf_id = [0]

names, loop_times, particles = get_data(ms,mpi,vf_id)

ave_time = []
run_name = []
body_count = []
p_s_count = []
time_per_particle = []
time_per_core = []
time_per_core_per_particle = []
time_per_core_per_bod = []

primary_sphere_count = [1,13,43,93,163,1,13,43,93,163]
core_count = [4,6,18,32]

for i, mpi_str in enumerate(mpi):
    ave_time.append([])
    run_name.append([])
    body_count.append([])
    p_s_count.append([])
    time_per_core.append([])
    time_per_particle.append([])
    for (n, stuf, numb, psc) in zip(names,loop_times[i], particles[i],primary_sphere_count):
        
        stufd = np.array(stuf)
        #print(np.mean(stuf))
        run_name[i].append(mpi_str + n)
        ave_time[i].append(np.mean(stufd[2:]))
        body_count[i].append(numb[0])
        p_s_count[i].append(numb[0] * psc)
        time_per_particle[i].append(np.mean(stufd[2:]) / (numb[0] * psc))
        time_per_core[i].append(np.mean(stufd[2:]) * (core_count[i]))
        time_per_core_per_particle.append(np.mean(stufd[2:]) / (core_count[i] * numb[0] * psc))
        time_per_core_per_bod.append(np.mean(stufd[2:]) / (core_count[i] * numb[0]))
   


x = np.linspace(np.min(p_s_count[i]), np.max(p_s_count[i]))
y = 0.00000001 * x**2
nlogn = 0.00000001 * x * np.log(x)

plt.figure(1)

plt.loglog(x,y,label="x^2", color='red', linestyle='dashed')
plt.loglog(x,nlogn,label="nlogn", color='black', linestyle='dashed')
for i, mpi_str in enumerate(mpi):
    plt.loglog(p_s_count[i],ave_time[i], marker="." , markersize= 12, linestyle="--", label=mpi_str + "vf0.05")


plt.xlabel(" Primary Spheres #")
plt.ylabel(" Loop Time (s) ")
plt.legend()
plt.savefig("low_vf_pcount_vs_time.png")

# plt.show()



temp_core_time = time_per_core
temp_psc = p_s_count

same_runs_tims_low = np.moveaxis(ave_time, [0,1], [1,0])
#same_runs_parts = np.moveaxis(p_s_count, [0,1], [1,0])
lables_low = ["s_0.05","m13_0.05","m43_0.05","m93_0.05","m163_0.05"]





ms = ["13_","43_","93_","163_"]
mpi = ["MPI221_","MPI321_","MPI332_","MPI442_"]
vf_id = [1]

names, loop_times, particles = get_data(ms,mpi,vf_id)

ave_time = []
run_name = []
body_count = []
p_s_count = []
time_per_particle = []
time_per_core = []
time_per_core_per_particle = []
time_per_core_per_bod = []

primary_sphere_count = [1,13,43,93,163,1,13,43,93,163]
core_count = [4,6,18,32]

for i, mpi_str in enumerate(mpi):
    ave_time.append([])
    run_name.append([])
    body_count.append([])
    p_s_count.append([])
    time_per_core.append([])
    time_per_particle.append([])
    for (n, stuf, numb, psc) in zip(names,loop_times[i], particles[i],primary_sphere_count):
        
        stufd = np.array(stuf)
        #print(np.mean(stuf))
        run_name[i].append(mpi_str + n)
        ave_time[i].append(np.mean(stufd[2:]))
        body_count[i].append(numb[0])
        p_s_count[i].append(numb[0] * psc)
        time_per_particle[i].append(np.mean(stufd[2:]) / (numb[0] * psc))
        time_per_core[i].append(np.mean(stufd[2:]) * (core_count[i]))
        time_per_core_per_particle.append(np.mean(stufd[2:]) / (core_count[i] * numb[0] * psc))
        time_per_core_per_bod.append(np.mean(stufd[2:]) / (core_count[i] * numb[0]))
   


x = np.linspace(np.min(p_s_count[i]), np.max(p_s_count[i]))
y = 0.00000001 * x**2
nlogn = 0.00000001 * x * np.log(x)

plt.figure(3)

plt.loglog(x,y,label="x^2", color='red', linestyle='dashed')
plt.loglog(x,nlogn,label="nlogn", color='black', linestyle='dashed')
for i, mpi_str in enumerate(mpi):
    plt.loglog(p_s_count[i],ave_time[i], marker="." , markersize= 12, linestyle="--", label=mpi_str + "vf0.5")


plt.xlabel(" Primary Spheres #")
plt.ylabel(" Loop Time (s) ")
plt.legend()
plt.savefig("high_vf_pcount_vs_time.png")
# plt.show()


plt.figure(2)

for i, mpi_str in enumerate(mpi):
    plt.loglog(p_s_count[i],time_per_core[i], marker="." , markersize= 12, linestyle="--",  label=mpi_str + "vf0.5")


for i, mpi_str in enumerate(mpi):
    plt.loglog(temp_psc[i],temp_core_time[i], marker="." , markersize= 12, linestyle="--",  label=mpi_str + "vf0.05")


plt.xlabel(" Primary Spheres #")
plt.ylabel(" Core Count x Loop Time (s) ")
plt.legend()
plt.savefig("core_time.png")
# plt.show()



same_runs_tims_high = np.moveaxis(ave_time, [0,1], [1,0])
same_runs_tims_core_high = np.moveaxis(time_per_core, [0,1], [1,0])
same_runs_tims_core_low = np.moveaxis(temp_core_time, [0,1], [1,0])
#same_runs_parts = np.moveaxis(p_s_count, [0,1], [1,0])
lables_high = ["s_0.5","m13_0.5","m43_0.5","m93_0.5","m163_0.5"]


plt.figure(25)
for i,same in enumerate(same_runs_tims_low):
    plt.semilogy(core_count, same, marker="." , markersize= 12, linestyle="--", label=lables_low[i])
plt.xlabel(" core #")
plt.ylabel(" Loop Time (s) ")
plt.xticks(core_count)
plt.legend()
plt.savefig("core_count_vs_time_low.png")

# plt.show()    

plt.figure(26)
for i,same in enumerate(same_runs_tims_high):
    if i < 3:
        if len(same) == 4:
            plt.semilogy(core_count, same, marker="." , markersize= 12, linestyle="--", label=lables_high[i])
        if len(same) == 3:
            plt.semilogy(core_count[:-1], same, marker="." , markersize= 12, linestyle="--", label=lables_high[i])
    # if len(same) == 2:
    #     plt.loglog(core_count[:-2], same, label=lables_high[i])
    # if len(same) == 1:
    #     plt.loglog(core_count[:-3], same, label=lables_high[i])
plt.xlabel(" core #")
plt.ylabel(" Loop Time (s) ")
plt.xticks(core_count)
plt.legend()
plt.savefig("core_count_vs_time_high.png")

# plt.show()    





plt.figure(27)
for i,same in enumerate(same_runs_tims_core_low):
    plt.semilogy(core_count, same, marker="." , markersize= 12, linestyle="--", label=lables_low[i])
plt.xlabel(" core #")
plt.ylabel("core # x Loop Time (s) ")
plt.xticks(core_count)
plt.legend()
plt.savefig("core_count_vs_timecore_low.png")

# plt.show()    

plt.figure(28)
for i,same in enumerate(same_runs_tims_core_high):
    if i < 3:
        if len(same) == 4:
            plt.semilogy(core_count, same, marker="." , markersize= 12, linestyle="--", label=lables_high[i])
        if len(same) == 3:
            plt.semilogy(core_count[:-1], same, marker="." , markersize= 12, linestyle="--", label=lables_high[i])
    # if len(same) == 2:
    #     plt.loglog(core_count[:-2], same, label=lables_high[i])
    # if len(same) == 1:
    #     plt.loglog(core_count[:-3], same, label=lables_high[i])
plt.xlabel(" core #")
plt.ylabel("core # x Loop Time (s) ")
plt.xticks(core_count)
plt.legend()
plt.savefig("core_count_vs_timecore_high.png")

# plt.show()    



