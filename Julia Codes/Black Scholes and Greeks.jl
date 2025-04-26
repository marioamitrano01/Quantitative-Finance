import Pkg
Pkg.add.(["Plots", "Distributions"])

using Plots
using Distributions
using Dates

function bs_d1(S, K, r, σ, T)
    return (log(S/K) + (r + 0.5*σ^2)*T) / (σ*sqrt(T))
end


function bs_d2(S, K, r, σ, T)
    return bs_d1(S, K, r, σ, T) - σ*sqrt(T)
end


function bs_call(S, K, r, σ, T)
    d1 = bs_d1(S, K, r, σ, T)
    d2 = bs_d2(S, K, r, σ, T)
    
    return S * cdf(Normal(), d1) - K * exp(-r*T) * cdf(Normal(), d2)
end


function bs_put(S, K, r, σ, T)
    d1 = bs_d1(S, K, r, σ, T)
    d2 = bs_d2(S, K, r, σ, T)
    
    return K * exp(-r*T) * cdf(Normal(), -d2) - S * cdf(Normal(), -d1)
end


function option_delta(S, K, r, σ, T, option_type::Symbol)
    d1 = bs_d1(S, K, r, σ, T)
    
    if option_type == :call
        return cdf(Normal(), d1)
    else # put
        return cdf(Normal(), d1) - 1.0
    end
end


function option_gamma(S, K, r, σ, T)
    d1 = bs_d1(S, K, r, σ, T)
    
    return pdf(Normal(), d1) / (S * σ * sqrt(T))
end


function option_vega(S, K, r, σ, T)
    d1 = bs_d1(S, K, r, σ, T)
    
    return 0.01 * S * sqrt(T) * pdf(Normal(), d1)
end


function option_theta(S, K, r, σ, T, option_type::Symbol)
    d1 = bs_d1(S, K, r, σ, T)
    d2 = bs_d2(S, K, r, σ, T)
    
    if option_type == :call
        theta_annual = -S * pdf(Normal(), d1) * σ / (2 * sqrt(T)) - 
                       r * K * exp(-r * T) * cdf(Normal(), d2)
    else 
        theta_annual = -S * pdf(Normal(), d1) * σ / (2 * sqrt(T)) + 
                       r * K * exp(-r * T) * cdf(Normal(), -d2)
    end
    
    return theta_annual / 365.0
end

function option_rho(S, K, r, σ, T, option_type::Symbol)
    d2 = bs_d2(S, K, r, σ, T)
    
    if option_type == :call
        return 0.01 * K * T * exp(-r * T) * cdf(Normal(), d2)
    else # put
        return -0.01 * K * T * exp(-r * T) * cdf(Normal(), -d2)
    end
end

function price_surface(option_type::Symbol, K, r, σ)
    S_range = range(0.5*K, 1.5*K, length=50)
    T_range = range(0.1, 1.0, length=50)
    
    prices = zeros(length(S_range), length(T_range))
    
    for (i, S) in enumerate(S_range)
        for (j, T) in enumerate(T_range)
            if option_type == :call
                prices[i, j] = bs_call(S, K, r, σ, T)
            else # put
                prices[i, j] = bs_put(S, K, r, σ, T)
            end
        end
    end
    
    return S_range, T_range, prices
end


function generate_price_surface_plot(option_type::Symbol, K, r, σ, output_dir)
    S_range, T_range, prices = price_surface(option_type, K, r, σ)
    
    option_name = string(option_type)
    plt = surface(S_range, T_range, prices, 
                 xlabel="Stock Price (S)", ylabel="Time to Maturity (T)", zlabel="Option Price",
                 title="$(uppercase(option_name)) Option Price Surface",
                 camera=(30, 30), color=:viridis,
                 size=(900, 700))
    
    savefig(plt, joinpath(output_dir, "$(option_name)_price_surface.png"))
    return plt
end


function generate_individual_greeks_plots(option_type::Symbol, K, r, σ, T, output_dir)
    S_range = range(0.5*K, 1.5*K, length=100)
    option_name = string(option_type)
    
    delta_values = [option_delta(S, K, r, σ, T, option_type) for S in S_range]
    delta_plot = plot(S_range, delta_values, 
                     xlabel="Stock Price", ylabel="Δ", 
                     title="$(uppercase(option_name)) Option Delta",
                     legend=false, linewidth=3, 
                     grid=true, framestyle=:box,
                     size=(800, 600))
    savefig(delta_plot, joinpath(output_dir, "$(option_name)_delta.png"))
    
    
    gamma_values = [option_gamma(S, K, r, σ, T) for S in S_range]
    gamma_plot = plot(S_range, gamma_values, 
                     xlabel="Stock Price", ylabel="Γ", 
                     title="$(uppercase(option_name)) Option Gamma",
                     legend=false, linewidth=3, 
                     grid=true, framestyle=:box,
                     size=(800, 600))
    savefig(gamma_plot, joinpath(output_dir, "$(option_name)_gamma.png"))
    
    vega_values = [option_vega(S, K, r, σ, T) / 0.01 for S in S_range]  
    vega_plot = plot(S_range, vega_values, 
                    xlabel="Stock Price", ylabel="ν", 
                    title="$(uppercase(option_name)) Option Vega (per 1% change)",
                    legend=false, linewidth=3, 
                    grid=true, framestyle=:box,
                    size=(800, 600))
    savefig(vega_plot, joinpath(output_dir, "$(option_name)_vega.png"))
    
    theta_values = [option_theta(S, K, r, σ, T, option_type) * 365 for S in S_range]  
    theta_plot = plot(S_range, theta_values, 
                     xlabel="Stock Price", ylabel="Θ", 
                     title="$(uppercase(option_name)) Option Theta (annual)",
                     legend=false, linewidth=3, 
                     grid=true, framestyle=:box,
                     size=(800, 600))
    savefig(theta_plot, joinpath(output_dir, "$(option_name)_theta.png"))
    
    
    rho_values = [option_rho(S, K, r, σ, T, option_type) / 0.01 for S in S_range]  
    rho_plot = plot(S_range, rho_values, 
                   xlabel="Stock Price", ylabel="ρ", 
                   title="$(uppercase(option_name)) Option Rho (per 1% change)",
                   legend=false, linewidth=3, 
                   grid=true, framestyle=:box,
                   size=(800, 600))
    savefig(rho_plot, joinpath(output_dir, "$(option_name)_rho.png"))
    
    return [delta_plot, gamma_plot, vega_plot, theta_plot, rho_plot]
end


function generate_volatility_plot(option_type::Symbol, S, K, r, T, output_dir)
    σ_range = range(0.05, 0.5, length=100)
    
    prices = zeros(length(σ_range))
    for (i, σ) in enumerate(σ_range)
        if option_type == :call
            prices[i] = bs_call(S, K, r, σ, T)
        else # put
            prices[i] = bs_put(S, K, r, σ, T)
        end
    end
    
    option_name = string(option_type)
    plt = plot(σ_range, prices, 
               xlabel="Volatility (σ)", ylabel="Option Price", 
               title="$(uppercase(option_name)) Option Price vs Volatility",
               legend=false, linewidth=3, 
               grid=true, framestyle=:box,
               size=(800, 600))
    
    savefig(plt, joinpath(output_dir, "$(option_name)_vs_volatility.png"))
    return plt
end


function generate_time_decay_plot(option_type::Symbol, S, K, r, σ, output_dir)
    T_range = range(0.01, 2.0, length=100)
    
    prices = zeros(length(T_range))
    for (i, T) in enumerate(T_range)
        if option_type == :call
            prices[i] = bs_call(S, K, r, σ, T)
        else # put
            prices[i] = bs_put(S, K, r, σ, T)
        end
    end
    
    option_name = string(option_type)
    plt = plot(T_range, prices, 
               xlabel="Time to Maturity (years)", ylabel="Option Price", 
               title="$(uppercase(option_name)) Option Price vs Time to Maturity",
               legend=false, linewidth=3, 
               grid=true, framestyle=:box,
               size=(800, 600))
    
    savefig(plt, joinpath(output_dir, "$(option_name)_time_decay.png"))
    return plt
end


function verify_formulas(S, K, r, σ, T)
    results = Dict()
    
    call_price = bs_call(S, K, r, σ, T)
    put_price = bs_put(S, K, r, σ, T)
    
    parity_lhs = call_price - put_price
    parity_rhs = S - K * exp(-r * T)
    parity_diff = abs(parity_lhs - parity_rhs)
    
    results["Call Price"] = call_price
    results["Put Price"] = put_price
    results["Put-Call Parity LHS"] = parity_lhs
    results["Put-Call Parity RHS"] = parity_rhs
    results["Put-Call Parity Difference"] = parity_diff
    
    call_delta = option_delta(S, K, r, σ, T, :call)
    put_delta = option_delta(S, K, r, σ, T, :put)
    delta_sum = call_delta + abs(put_delta)
    
    results["Call Delta"] = call_delta
    results["Put Delta"] = put_delta
    results["Delta Relationship"] = delta_sum  
    
    call_gamma = option_gamma(S, K, r, σ, T)
    put_gamma = option_gamma(S, K, r, σ, T)
    gamma_diff = abs(call_gamma - put_gamma)
    
    call_vega = option_vega(S, K, r, σ, T)
    put_vega = option_vega(S, K, r, σ, T)
    vega_diff = abs(call_vega - put_vega)
    
    results["Call Gamma"] = call_gamma
    results["Put Gamma"] = put_gamma
    results["Gamma Difference"] = gamma_diff  
    
    results["Call Vega"] = call_vega
    results["Put Vega"] = put_vega
    results["Vega Difference"] = vega_diff 
    
    return results
end


function execute_black_scholes_analysis(S, K, r, σ, T)
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    output_dir = joinpath(homedir(), "Desktop", "BlackScholes_Analysis_")
    mkpath(output_dir)
    
    verification = verify_formulas(S, K, r, σ, T)
    
    open(joinpath(output_dir, "analysis_summary.txt"), "w") do io
        write(io, "BLACK-SCHOLES ANALYSIS SUMMARY\n")
        write(io, "PARAMETERS\n")
        write(io, "S = $(S)\n")
        write(io, "K = $(K)\n")
        write(io, "r = $(r)\n")
        write(io, "σ = $(σ)\n")
        write(io, "T = $(T)\n\n")
        
        write(io, "RESULTS\n")
        write(io, "Call Price: $(verification["Call Price"])\n")
        write(io, "Put Price: $(verification["Put Price"])\n")
        
        write(io, "\nGREEKS\n")
        write(io, "Call Delta: $(verification["Call Delta"])\n")
        write(io, "Put Delta: $(verification["Put Delta"])\n")
        write(io, "Gamma: $(verification["Call Gamma"])\n")
        write(io, "Vega (per 1%): $(verification["Call Vega"]/0.01)\n")
        write(io, "Call Theta (annual): $(option_theta(S, K, r, σ, T, :call) * 365)\n")
        write(io, "Put Theta (annual): $(option_theta(S, K, r, σ, T, :put) * 365)\n")
        write(io, "Call Rho (per 1%): $(option_rho(S, K, r, σ, T, :call)/0.01)\n")
        write(io, "Put Rho (per 1%): $(option_rho(S, K, r, σ, T, :put)/0.01)\n")
        
        write(io, "\nVERIFICATION\n")
        write(io, "Put-Call Parity: C - P = S - Ke^(-rT)\n")
        write(io, "                 $(verification["Put-Call Parity LHS"]) = $(verification["Put-Call Parity RHS"])\n")
        write(io, "                 Difference: $(verification["Put-Call Parity Difference"])\n\n")
        
        write(io, "Delta Relationship: Δ_call + |Δ_put| ≈ 1\n")
        write(io, "                    $(verification["Call Delta"]) + |$(verification["Put Delta"])| = $(verification["Delta Relationship"])\n\n")
        
        write(io, "Gamma Equality: Γ_call = Γ_put\n")
        write(io, "                $(verification["Call Gamma"]) = $(verification["Put Gamma"])\n")
        write(io, "                Difference: $(verification["Gamma Difference"])\n\n")
        
        write(io, "Vega Equality: ν_call = ν_put\n")
        write(io, "              $(verification["Call Vega"]) = $(verification["Put Vega"])\n")
        write(io, "              Difference: $(verification["Vega Difference"])\n")
    end
    
    generate_price_surface_plot(:call, K, r, σ, output_dir)
    generate_price_surface_plot(:put, K, r, σ, output_dir)
    generate_individual_greeks_plots(:call, K, r, σ, T, output_dir)
    generate_individual_greeks_plots(:put, K, r, σ, T, output_dir)
    generate_volatility_plot(:call, S, K, r, T, output_dir)
    generate_volatility_plot(:put, S, K, r, T, output_dir)
    generate_time_decay_plot(:call, S, K, r, σ, output_dir)
    generate_time_decay_plot(:put, S, K, r, σ, output_dir)
    
    open(joinpath(output_dir, "formulas.txt"), "w") do io
        write(io, "BLACK-SCHOLES MODEL FORMULAS\n")
        write(io, "CORE EQUATIONS\n")
        write(io, "d₁ = (ln(S/K) + (r + σ²/2)T) / (σ√T)\n")
        write(io, "d₂ = d₁ - σ√T\n\n")
        
        write(io, "OPTION PRICES\n")
        write(io, "Call Price = S·N(d₁) - Ke⁻ʳᵀ·N(d₂)\n")
        write(io, "Put Price = Ke⁻ʳᵀ·N(-d₂) - S·N(-d₁)\n\n")
        
        write(io, "WHERE\n")
        write(io, "S = Current stock price\n")
        write(io, "K = Strike price\n")
        write(io, "r = Risk-free interest rate (annual, continuously compounded)\n")
        write(io, "σ = Volatility of the stock (annual)\n")
        write(io, "T = Time to expiration (in years)\n")
        write(io, "N(x) = Cumulative distribution function of the standard normal distribution\n\n")
        
        write(io, "GREEKS\n")
        write(io, "Call Delta = N(d₁)\n")
        write(io, "Put Delta = N(d₁) - 1\n")
        write(io, "Gamma = N'(d₁) / (S·σ·√T)  (same for calls and puts)\n")
        write(io, "Vega = S·√T·N'(d₁)  (same for calls and puts)\n")
        write(io, "Call Theta = -S·N'(d₁)·σ/(2·√T) - r·K·e⁻ʳᵀ·N(d₂)\n")
        write(io, "Put Theta = -S·N'(d₁)·σ/(2·√T) + r·K·e⁻ʳᵀ·N(-d₂)\n")
        write(io, "Call Rho = K·T·e⁻ʳᵀ·N(d₂)\n")
        write(io, "Put Rho = -K·T·e⁻ʳᵀ·N(-d₂)\n\n")
        
        write(io, "KEY RELATIONSHIPS\n")
        write(io, "Put-Call Parity: C - P = S - Ke⁻ʳᵀ\n")
        write(io, "Delta Relationship: Δ_call + |Δ_put| = 1\n")
        write(io, "Gamma Equality: Γ_call = Γ_put\n")
        write(io, "Vega Equality: ν_call = ν_put\n")
    end
    
    return output_dir
end

S = 100.0
K = 100.0
r = 0.05
σ = 0.2
T = 1.0

output_dir = execute_black_scholes_analysis(S, K, r, σ, T)
println("Analysis complete. Results saved in: $(output_dir)")
