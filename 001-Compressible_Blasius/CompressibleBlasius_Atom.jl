using BenchmarkTools
using Plots
using Printf
using DelimitedFiles
using LinearAlgebra

function Y1(y₂)
    return y₂
end

function Y2(y₃)
    return y₃
end

function Y3(y₁,y₃,y₄,y₅,cμ,T∞)
    return -y₃*((y₅/(2*(y₄)))-(y₅/(y₄+cμ/T∞)))-y₁*y₃*((y₄+cμ/T∞)/(sqrt(y₄)*(1+cμ/T∞)));
end

function Y4(y₅)
    return y₅
end

function Y5(y₁,y₃,y₄,y₅,cμ,T∞,M∞,Pr,γ)
    return -y₅^2*((0.5/y₄)-(1/(y₄+cμ/T∞)))-Pr*y₁*y₅/sqrt(y₄)*(y₄+cμ/T∞)/(1+cμ/T∞)-(γ-1)*Pr*M∞^2*y₃^2;
end

function RK(N,h,y1,y2,y3,y4,y5,cμ,T∞,Pr,γ,M∞)
	for i=1:N

		k11 = Y1(y2[i]);
		k21 = Y2(y3[i]);
		k31 = Y3(y1[i], y3[i], y4[i], y5[i],cμ,T∞);
		k41 = Y4(y5[i]);
		k51 = Y5(y1[i], y3[i], y4[i], y5[i],cμ,T∞,M∞,Pr,γ);

		k12 = Y1(y2[i]+0.5*h*k21);
		k22 = Y2(y3[i]+0.5*h*k31);
		k32 = Y3(y1[i]+0.5*h*k11, y3[i]+0.5*h*k31, y4[i]+0.5*h*k41, y5[i]+0.5*h*k51,cμ,T∞);
		k42 = Y4(y5[i]+0.5*h*k51);
		k52 = Y5(y1[i]+0.5*h*k11, y3[i]+0.5*h*k31, y4[i]+0.5*h*k41, y5[i]+0.5*h*k51,cμ,T∞,M∞,Pr,γ);

		k13 = Y1(y2[i]+0.5*h*k22);
		k23 = Y2(y3[i]+0.5*h*k32);
		k33 = Y3(y1[i]+0.5*h*k12, y3[i]+0.5*h*k32, y4[i]+0.5*h*k42, y5[i]+0.5*h*k52,cμ,T∞);
		k43 = Y4(y5[i]+0.5*h*k52);
		k53 = Y5(y1[i]+0.5*h*k12, y3[i]+0.5*h*k32, y4[i]+0.5*h*k42, y5[i]+0.5*h*k52,cμ,T∞,M∞,Pr,γ);

		k14 = Y1(y2[i]+h*k23);
		k24 = Y2(y3[i]+h*k33);
		k34 = Y3(y1[i]+h*k13, y3[i]+h*k33, y4[i]+h*k43, y5[i]+h*k53,cμ,T∞);
		k44 = Y4(y5[i]+h*k53);
		k54 = Y5(y1[i]+h*k13, y3[i]+h*k33, y4[i]+h*k43, y5[i]+h*k53,cμ,T∞,M∞,Pr,γ);

		y5[i+1] = y5[i] + (1/6)*(k51 + 2*k52 + 2*k53 + k54)*h;
		y4[i+1] = y4[i] + (1/6)*(k41 + 2*k42 + 2*k43 + k44)*h;
		y3[i+1] = y3[i] + (1/6)*(k31 + 2*k32 + 2*k33 + k34)*h;
		y2[i+1] = y2[i] + (1/6)*(k21 + 2*k22 + 2*k23 + k24)*h;
		y1[i+1] = y1[i] + (1/6)*(k11 + 2*k12 + 2*k13 + k14)*h;
	end
	return y1,y2,y3,y4,y5
end

function selfsimilar(M∞=1.0, T∞=300, Tw=2, ηmax=10, N=101, ϵ=1e-9)
    Δη = ηmax/N
    γ  = 1.4   # Ratio of specific heat
    cμ = 110.4 # Sutherland law coefficient for [Kelvin]
    Pr = 0.72  # Prandtl Number
	adi = 1

    Δ = 1e-10     # Small number for shooting method (Decrease when ϵ is decreased)

    # Initializing the solution vectors
    y1 = zeros(N+1)   # f
    y2 = zeros(N+1)   # f'
    y3 = zeros(N+1)   # f''
    y4 = zeros(N+1)   # ρ(η)
    y5 = zeros(N+1)   # ρ(η)'
    η  = [i*Δη for i=0:N]

	if adi==1
		y1[1] = 0.0
    	y2[1] = 0.0
    	y5[1] = 0.0

    	α₀ = 0.1 	  # Initial Guess
    	β₀ = 2.3118      # Initial Guess
	elseif adi==0
		y1[1] = 0
    	y2[1] = 0
    	y4[1] = Tw

    	α₀ = 0.1 	  # Initial Guess
    	β₀ = 3.0      # Initial Guess
	end

    y2n = 1
    y4n = 1
    for i=1:100000
        if adi==1
    		y1[1] = 0.0
        	y2[1] = 0.0
        	y5[1] = 0.0

            y3[1] = α₀;    # Initial Guess
            y4[1] = β₀;    # Initial Guess
    	elseif adi==0
    		y1[1] = 0
        	y2[1] = 0
        	y4[1] = Tw

            y3[1] = α₀;    # Initial Guess
            y5[1] = β₀;    # Initial Guess
    	end


		# First solution for Newton's iteration
        y1,y2,y3,y4,y5 = RK(N,Δη,y1,y2,y3,y4,y5,cμ,T∞,Pr,γ,M∞)

        # Storing the freestream values for Newton's iteration method
        y₂ₒ = y2[N+1];
        y₄ₒ = y4[N+1];

        # Small number addition for Newton's iteration method
        if adi==1
    		y1[1] = 0.0
        	y2[1] = 0.0
        	y5[1] = 0.0

            y3[1] = α₀+Δ;   # Initial Guess + Small number
            y4[1] = β₀;    # Initial Guess
    	elseif adi==0
    		y1[1] = 0
        	y2[1] = 0
        	y4[1] = Tw

            y3[1] = α₀+Δ;    # Initial Guess + Small number
            y5[1] = β₀;    # Initial Guess
    	end

		# Second solution for Newton's iteration
        y1,y2,y3,y4,y5 = RK(N,Δη,y1,y2,y3,y4,y5,cμ,T∞,Pr,γ,M∞)

        # Storing the freestream values for Newton's iteration method
        y₂ₙ₁ = y2[N+1];
        y₄ₙ₁ = y4[N+1];

        # Small number addition for Newton's iteration method
        # Small number addition for Newton's iteration method
        if adi==1
    		y1[1] = 0.0
        	y2[1] = 0.0
        	y5[1] = 0.0

            y3[1] = α₀;   # Initial Guess
            y4[1] = β₀+Δ;    # Initial Guess + Small number
    	elseif adi==0
    		y1[1] = 0
        	y2[1] = 0
        	y4[1] = Tw

            y3[1] = α₀;    # Initial Guess
            y5[1] = β₀+Δ;    # Initial Guess + Small number
    	end

		# Third solution for Newton's iteration
        y1,y2,y3,y4,y5 = RK(N,Δη,y1,y2,y3,y4,y5,cμ,T∞,Pr,γ,M∞)

        # Storing the freestream values for Newton's iteration method
        y₂ₙ₂ = y2[N+1];
        y₄ₙ₂ = y4[N+1];

        # Calculation of the next initial guess with Newton's iteration method
        p₁₁ = (y₂ₙ₁-y₂ₒ)/Δ;
        p₂₁ = (y₄ₙ₁-y₄ₒ)/Δ;
        p₁₂ = (y₂ₙ₂-y₂ₒ)/Δ;
        p₂₂ = (y₄ₙ₂-y₄ₒ)/Δ;
        r₁ = 1-y₂ₒ;
        r₂ = 1-y₄ₒ;
        Δα = (p₂₂*r₁-p₁₂*r₂)/(p₁₁*p₂₂-p₁₂*p₂₁);
        Δβ = (p₁₁*r₂-p₂₁*r₁)/(p₁₁*p₂₂-p₁₂*p₂₁);
        α₀ = α₀ + Δα;
        β₀ = β₀ + Δβ;

        #@show(norm(y4))
		if (abs(y2[N+1]-1.0)<ϵ) && (abs(y4[N+1]-1.0)<ϵ) && (abs(y2n-norm(y2))<ϵ) && (abs(y4n-norm(y4))<ϵ)
        	break
    	end
        y2n = norm(y2)
        y4n = norm(y4)
    end

	# Copying values for logical names
    U = y2
    T = y4

    # Integration for η --> y transformation
    y = zeros(N+1);
    for i=2:N+1
       y[i] = y[i-1] + y4[i]*(η[i]-η[i-1]);
    end
    y = y*sqrt(2);

    return η,y,U,T,N
end

# ╔═╡ f044ac30-60d8-4ebe-8048-463ba4b8bebd
ηtest, ytest, Utest, Ttest, Ntest = selfsimilar();
result = @benchmark selfsimilar()

plot([Utest,Ttest],[ηtest,ηtest],
        title = "Compressible Blasius Profiles",
        label = ["U" "T"],
        legend = :topleft,
        xlabel = "U,T",
        ylabel = "\\eta",
        linewidth = 2,
        linecolor = :black,
        markershape = :circle,
        markercolor = :auto,
	)
