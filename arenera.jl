using Random

function getTime(prob,attempts_per_t,num_samples)
    t_arr = zeros(num_samples)
    for i in 1:num_samples
        t = 1
        while all([rand() > prob for _ in 1:attempts_per_t])
            t +=1
        end
        t_arr[i] = t
    end

    return t_arr
end

function getTimeSplitRate(prob,attempts_per_t,num_samples,split,split_t)
    t_change = split_t * 1.0/prob
    println((t_change,prob*split[1]/(split_t*sum(split)),prob*split[2]/(split_t*sum(split))))
    t_arr = zeros(num_samples)
    for i in 1:num_samples
        t = 1
        while all([rand() > prob*split[1]/(split_t*sum(split)) for _ in 1:attempts_per_t]) && t < t_change
            t +=1
        end
        while all([rand() > prob*split[2]/(split_t*sum(split)) for _ in 1:attempts_per_t]) && t >= t_change
            t +=1
        end
        t_arr[i] = t
    end

    return t_arr
end

function getTimeSplitAttempts(prob,attempts_per_t,num_samples,split,split_t)
    # att_1 = attempts_per_t*split[1]*2/sum(split)
    # att_2 = attempts_per_t*split[2]*2/sum(split)
    # att_1 = att_1 / min(att_1,att_2)
    # att_2 = floor(att_2 / min(att_1,att_2))
    t_change = split_t * 1.0/prob
    println((t_change,split[1]/(split_t*sum(split)),split[2]/(split_t*sum(split))))
    t_arr = zeros(num_samples)
    for i in 1:num_samples
        t = 1
        while all([rand() > prob for _ in 1:attempts_per_t]) && t < t_change
            t +=1
        end
        t1 = t / (split[1]/(split_t*sum(split)))
        t = 0
        while all([rand() > prob for _ in 1:attempts_per_t]) && t >= t_change
            t +=1
        end
        t2 = t / (split[2]/(split_t*sum(split)))
        t_arr[i] = t1 + t2
    end

    return t_arr
end

n=100000
p=0.01
f=10
xmax=1000

ts1 = getTime(p,1,n);
histogram(ts1,xlimits=(0,xmax))

ts2 = getTime(p,f,n);
histogram(ts2,xlimits=(0,xmax))

ts3 = getTime(f*p,1,n);
histogram(ts3,xlimits=(0,xmax))

println((mean(ts1),mean(ts2),mean(ts3),mean(ts2)/mean(ts1),mean(ts3)/mean(ts1)))

first_t = minimum(hcat(ts1,ts2,ts3),dims=2);
println((mean(first_t),1.0/(1.0/mean(ts1)+1.0/mean(ts2)+1.0/mean(ts3))))
evt_1 = (first_t.==ts1);
evt_2 = (first_t.==ts2);
evt_3 = (first_t.==ts3);
println((mean(evt_1),mean(evt_2),mean(evt_3)))


ts4 = getTimeSplitRate(p,1,n,[2,1],0.5);
histogram(ts4,xlimits=(0,xmax))

println((mean(ts1),mean(ts2),mean(ts3),mean(ts4),mean(ts2)/mean(ts1),mean(ts3)/mean(ts1),mean(ts4)/mean(ts1)))

first_t = minimum(hcat(ts1,ts4),dims=2);
println((mean(first_t),1.0/(1.0/mean(ts1)+1.0/mean(ts4))))
evt_1 = (first_t.==ts1);
evt_4 = (first_t.==ts4);
println((mean(evt_1),mean(evt_4)))
println((mean(evt_1 .* first_t),mean(evt_4 .* first_t)))

ts4 = getTimeSplitAttempts(p,1,n,[1,1],0.5);
histogram(ts4,xlimits=(0,xmax))

println((mean(ts1),mean(ts2),mean(ts3),mean(ts4),mean(ts2)/mean(ts1),mean(ts3)/mean(ts1),mean(ts4)/mean(ts1)))

first_t = minimum(hcat(ts1,ts4),dims=2);
println((mean(first_t),1.0/(1.0/mean(ts1)+1.0/mean(ts4))))
evt_1 = (first_t.==ts1);
evt_4 = (first_t.==ts4);
println((mean(evt_1),mean(evt_4)))




function gil(inf1,inf2,mut1,mut2,tf)
    t = 0
    vars = [0.0,0.0] # infected 1 and 2
    rates = [inf1-inf1*mut1,inf2-inf2*mut2,inf1*mut1,inf2*mut2] # infection 1 and 2, mutation 1 and 2
    r_tot = sum(rates)
    r_evt = r_tot * rand()
    t_next = (1.0/r_tot) * randexp()

    switched = false

    while t < tf
        rates = [inf1-inf1*mut1,inf2-inf2*mut2,inf1*mut1,inf2*mut2]
        r_tot = sum(rates)
        r_evt = r_tot * rand()
        t_next = (1.0/r_tot) * randexp()
        t += t_next
        if t > tf/2.0 && !switched
            t = tf/2.0
            inf1 = 2*inf1
            switched = true
        elseif r_evt < rates[1]
            vars[1] += 1
        elseif r_evt < sum(rates[1:2])
            vars[2] += 1
        elseif r_evt < sum(rates[1:3])
            return [1,t]
        elseif r_evt < sum(rates[1:4])
            return [2,t]
        else
            println("ERROR")
        end
    end

    return [0,t]
end

function gil2(inf1,inf2,mut1,mut2,tf)
    t = 0
    vars = [0.0,0.0] # infected 1 and 2
    rates = [inf1-inf1*mut1,inf2-inf2*mut2,inf1*mut1,inf2*mut2] # infection 1 and 2, mutation 1 and 2
    r_tot = sum(rates)
    r_evt = r_tot * rand()
    t_next = (1.0/r_tot) * randexp()

    switched = false

    while t < tf
        rates = [inf1-inf1*mut1,inf2-inf2*mut2,inf1*mut1,inf2*mut2]
        r_tot = sum(rates)
        r_evt = r_tot * rand()
        t_next = (1.0/r_tot) * randexp()
        t += t_next
        if t > tf/2.0 && !switched
            t = tf/2.0
            inf1 = 2*inf1
            switched = true
        elseif r_evt < rates[1]
            vars[1] += 1
        elseif r_evt < sum(rates[1:2])
            vars[2] += 1
        elseif r_evt < sum(rates[1:3])
            return [1,t]
        elseif r_evt < sum(rates[1:4])
            return [2,t]
        else
            println("ERROR")
        end
    end

    return [0,t]
end

gil(0.01,0.01,0.01,0.01,1000000)
