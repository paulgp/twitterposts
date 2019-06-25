

cd "~/Dropbox/TwitterPosts/"

set matsize 10000
foreach N in 400 1000 5000 {
	foreach b in 0.5 1.0 {
		clear all
		forvalues bs = 1/1 {
			clear
			disp "`bs'.."_continue
			quietly set obs `N'
			gen x = (runiform()*2) -1

        *Var(eps) = 1
			gen eps=rnormal()

         *True beta 
			gen treatment=`b'
			quietly replace treatment=0 if x < 0
			if "`b'" == "0.1" {
				local treatment_label "0_1"
				}
			else if "`b'" == "0.5" {
				local treatment_label "0_5"
				}
			else {
				local treatment_label "1_0"
				}
			gen outcome=treatment+ eps
			gen above= x>=0
			gen above_x = above*x
* Var(Y | x = 0) = Var(eps) = 1

*get p-value
			qui reg outcome above x above_x   , r
			mat table = r(table)
			local beta_wrong = table[1,1]
			local p_wrong = table[4,1]
			qui rdrobust outcome x
* estimated effect
* e(b)
			mat b = e(b)
			local beta = b[1,1]
			local p_val_robust = string(e(pv_rb), "%9.5f")
			local p_wrong = string(`p_wrong', "%9.5f")
* binned plot with a linear fit
			if `bs' == 100 | `bs' == 200 | `bs' == 300 {
				qui rdplot outcome x , p(1) ///
				  graph_options(subtitle("BetaHat = `beta', p-value = `p_val_robust', wrong p-value = `p_wrong'"))
				graph export "plot_`N'_`treatment_label'_`bs'.pdf", replace
				}
			scatter outcome x 
			qui count if inrange(x, -.2, .2)
			local num_obs = r(N)

			mat beta_mat = (nullmat(beta_mat), `beta')
			mat p_mat = (nullmat(p_mat), `p_val_robust')
			mat p_wrong_mat = (nullmat(p_wrong_mat), `p_wrong')
			mat obs_mat = (nullmat(obs_mat), `num_obs')		
			}


		clear

		foreach x in beta p obs p_wrong {
			mat `x'_mat = `x'_mat'
			svmat `x'_mat
			}

		hist beta_mat1, xtitle("Beta distribution") note("Sims = 400")
		graph export "plot_beta_`N'_`treatment_label'.pdf", replace
		hist p_mat1, xtitle("p-val distribution") note("Sims = 400") 
		graph export "plot_p_`N'_`treatment_label'.pdf", replace
		hist p_wrong_mat1, xtitle("Wrong p-val distribution") note("Sims = 400") 
		graph export "plot_pwrong_`N'_`treatment_label'.pdf", replace
		hist obs_mat1, xtitle("Effective observations distribution") note("Sims = 400") 
		graph export "plot_obs_`N'_`treatment_label'.pdf", replace
		}
	}
/*** output **/
/* Thre example graphs
* report beta, N, effective num obs
* a
* b
* c
* beta hist
* p_val hist

**/
