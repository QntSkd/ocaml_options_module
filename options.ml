module Options =
struct 
  let erf x =
    let a1 = 0.254829592 in
    let a2 = -0.2844967368 in
    let a3 = 1.421413741 in
    let a4 = -1.453152027 in
    let a5 = 1.061405429 in
    let p = 0.3275911 in
    let sign = x/.Float.abs(x) in
    let x_1 = Float.abs(x) in
    let t = 1.0/.(1.0 +. p*.x_1) in
    let y = 1.0 -. (((((a5*.t +. a4)*.t) +. a3)*.t +. a2)*.t +. a1)*.t*.Float.exp(-1.*.x*.x) in
    sign *. y
  
  let npdf x =
    1.0/.Float.sqrt(2.0*.Float.pi)*.Float.exp(-0.5*.(x*.x))
  
  let ncdf x =
    0.5 *. (1. +. erf (x/.Float.sqrt(2.)))
  
  let black_scholes s k t r v =
    let d1 = (Float.log(s/.k) +. (r +. (v*.v)/.2.)*.t)/.(v*.Float.sqrt(t)) in
    let d2 = d1 -. v *. Float.sqrt(t) in
    s *. ncdf d1 -. ncdf d2 *. k *. Float.exp(-.r *. t) 

  let vega s k t r v =
    let d1 = (Float.log(s/.k) +. (r +. (v*.v)/.2.0)*.t)/.v *.Float.sqrt(t) in
    s*. Float.sqrt(t) *. (npdf d1) 

  let rec implied_vol c s k t r v epsilon = 
    let diff = (black_scholes s k t r v) -. c in
    let veg = vega s k t r v in 
    match diff with 
    | diff when Float.abs(diff) < epsilon -> v
    | _ -> implied_vol c s k t r (v -. diff/.veg) epsilon
end
(* 
for black scholes call pricing:
    s - spot price
    k - strike price
    r - risk-free rate
    t - time until expiration
    v - volatility
    
for implied volatility calculation:
    s - spot price
    k - strike price
    r - risk-free rate
    t - time until expiration
    v - **implied volatility guess**
    epsilon - accuracy tolerance
*)

let s = 451.66
let k = 452.
let t = 0.05
let r = 0.016
let v = 0.2374
let call_price = Options.black_scholes s k t r v
let epsilon = 0.0001
let iv_calc = Options.implied_vol call_price s k t r v epsilon 
Printf.printf "Call Price: $%f\nKnown implied volatility: %f\nDerived implied volatility: %f\n" call_price v iv_calc
