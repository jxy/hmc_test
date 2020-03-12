using Random, SpecialFunctions, Statistics, UnicodePlots

struct OneLink
	beta
end
OneLink(;beta = 2.0) = OneLink(beta)
action(m::OneLink, x) = -m.beta*cos(x)
force(m::OneLink, x) = m.beta*sin(x)

struct HMC
	tau
	nstep
	dt
	seed
end
function HMC(;tau = 4.0, nstep = 10, seed = 11*13)
	dt = tau/nstep
	HMC(tau,nstep,dt,seed)
end

function leapfrog(hmc::HMC, m, x, p)
	dt = hmc.dt
	x_ = x + 0.5dt * p
	p_ = p + (-dt) * force(m, x_)
	for i = 1:hmc.nstep
		x_ += dt * p_
		p_ += (-dt) * force(m, x_)
	end
	x_ += 0.5dt * p_
    x_ = mod(x_, 2pi)
	(x_, -p_)	# flip p for reversibility
end

struct L2HMC
	tau
	nstep
	dt
	wv
	bv
	seed
end
function L2HMC(;tau = 4.0, nstep = 10, seed = 11*13, wv = sqrt(43), bv = sqrt(41))
	dt = tau/nstep
	L2HMC(tau,nstep,dt,wv,bv,seed)
end

tx(l2hmc::L2HMC, p) = p * l2hmc.wv + l2hmc.bv
function leapfrog(l2hmc::L2HMC, m, x, p, d)
	dt = l2hmc.dt
	if d > 0
		x_ = x + 0.5dt * (p + tx(l2hmc, p))
		p_ = p + (-dt) * force(m, x_)
		for i = 1:l2hmc.nstep
			x_ += dt * (p_ + tx(l2hmc, p_))
			p_ += (-dt) * force(m, x_)
		end
		x_ += 0.5dt * (p_ + tx(l2hmc, p_))
	else
		x_ = x - 0.5dt * (p + tx(l2hmc, p))
		p_ = p - (-dt) * force(m, x_)
		for i = 1:l2hmc.nstep
			x_ -= dt * (p_ + tx(l2hmc, p_))
			p_ -= (-dt) * force(m, x_)
		end
		x_ -= 0.5dt * (p_ + tx(l2hmc, p_))
	end
    x_ = mod(x_, 2pi)
	(x_, p_, -d)	# flip d for reversibility
end

struct HMCState
	x
	p
end
struct L2HMCState
	x
	p
	d
end

import Base.show
show(io::IO, st::HMCState) = print(io, "HMCState: $(st.x) $(st.p)")
show(io::IO, st::L2HMCState) = print(io, "L2HMCState: $(st.x) $(st.p) $(st.d)")

function fresh(mc::HMC, x)
	p = randn(eltype(x))
	HMCState(x, p)
end

function propose(hmc::HMC, m, st::HMCState)
    x, p = st.x, st.p
    x_, p_ = leapfrog(hmc, m, x, p)
	HMCState(x_, p_)
end

function fresh(mc::L2HMC, x)
	p = randn(eltype(x))
	d = rand() >= 0.5 ? 1 : -1
	L2HMCState(x, p, d)
end

function propose(l2hmc::L2HMC, m, st::L2HMCState)
    x, p, d = st.x, st.p, st.d
    x_, p_, d_ = leapfrog(l2hmc, m, x, p, d)
	L2HMCState(x_, p_, d_)
end

function gen(mc, m, st; io = nothing)
	act0 = action(m, st.x) + 0.5st.p*st.p
	pst = propose(mc, m, st)
	st_ = propose(mc, m, pst)
	dx = abs(st_.x-st.x)
	dp = abs(st_.p-st.p)
	if dx > 1e-10 || dp > 1e-10
		println("failed rev check: $st")
		println("	-> $pst")
		println("	<- $st_")
	end
    act = action(m, pst.x) + 0.5pst.p*pst.p
    prob = rand()
    dH = act-act0
    exp_mdH = exp(-dH)
    acc = prob < exp_mdH
    newst = acc ? pst : st
	if !isnothing(io)
		println(io, "$st $pst $newst")
	end
    (dx, dp, dH, exp_mdH, acc, newst.x)
end

genx(mc, m, x) = gen(mc, m, fresh(mc, x))

function mcmc(mc; ntraj = 1_000_000)
	Random.seed!(mc.seed)
	m = OneLink()
	cosx = besseli(1, m.beta)/besseli(0, m.beta)
	xinit = 3.0		# acos(cosx)
	xs = Array{Float64}(undef, ntraj)
	dxs = similar(xs)
	dps = similar(xs)
	dhs = similar(xs)
	exp_mdhs = similar(xs)
	accs = similar(xs)
	for i = 1:length(xs)
		(dxs[i], dps[i], dhs[i], exp_mdhs[i], accs[i], xs[i]) = genx(mc, m, i==1 ? xinit : xs[i-1])
	end
	println(mc)
	println(scatterplot(mod.(pi .+ xs, 2pi) .- pi, width=100, height=30, ylim=[-pi,pi]))
	# pick the later half for statistics
	nstati = ntraj + 1 - ntraj ÷ 2
	dxs = dxs[nstati:end]
	dps = dps[nstati:end]
	dhs = dhs[nstati:end]
	exp_mdhs = exp_mdhs[nstati:end]
	accs = accs[nstati:end]
	xs = xs[nstati:end]
	cosxs = cos.(xs)
	t = typeof(mc)
	nstat = length(cosxs)
	println("# using the last $nstat configurations")
	println("$t: <cos(x)-exact> ", mean(cosxs)-cosx, " +/- ", std(cosxs)/sqrt(nstat))
	println("$t: <dx> ", mean(dxs), " +/- ", std(dxs)/sqrt(nstat))
	println("$t: <dp> ", mean(dps), " +/- ", std(dps)/sqrt(nstat))
	println("$t: <dH> ", mean(dhs), " +/- ", std(dhs)/sqrt(nstat))
	#println("$t: <exp(-dH)> ", mean(exp_mdhs), " +/- ", std(exp_mdhs)/sqrt(nstat))
	println("$t: <acc> ", mean(accs), " +/- ", std(accs)/sqrt(nstat))
end

function test_mcmc()
	for i in [1,2,4,8,16,32]
		ntraj = i*1_000_000
		println("Ntraj = $ntraj")
		mcmc(HMC(tau=2.0, nstep=1, seed=13^7+i), ntraj = ntraj)
		mcmc(L2HMC(tau=2.0, nstep=1, wv = 0.1-1.0, bv = 3pi, seed=13^7+i), ntraj = ntraj)
		println()
	end
end

function test_ergo()
	nx = 200
	np = 200
	xs = range(pi/nx, 2pi, step = 2pi/nx)
	ps = range(-6, 6, length = np)
	m = OneLink()
	open("ergo_hmc.out", "w") do o
		mc = HMC(tau=2.0, nstep=1, seed=13^7)
		for x in xs
			for p in ps
				gen(mc, m, HMCState(x, p), io = o)
			end
		end
	end
	open("ergo_l2hmc.out", "w") do o
		mc = L2HMC(tau=2.0, nstep=1, wv = 0.1-1.0, bv = 3pi, seed=13^7)
		for d in [-1, 1]
			for x in xs
				for p in ps
					gen(mc, m, L2HMCState(x, p, d), io = o)
				end
			end
		end
	end
end

#test_mcmc()
test_ergo()
