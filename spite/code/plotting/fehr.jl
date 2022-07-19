
p=plot(xaxis = nothing, yaxis = nothing, right_margin = 20px, top_margin = 20px, left_margin = 60px)
xaxis!("xⱼ", draw_arrow = true, guidefonthalign = :right,foreground_color_border = :white)
yaxis!("", draw_arrow = true, guidefontvalign = :top, foreground_color_border = :white)
annotate!((-0.05,.9), "Uᵢ(xⱼ|xᵢ)",  9)

plot!([0,0],[0,1.1],arrow=true,color=:black,linewidth=2,label="")
plot!([0,1.1],[0,0],arrow=true,color=:black,linewidth=2,label="")
plot!([.5,.1],[.5,.4],color=:black,linewidth=1,label=false)
plot!([.5,.8],[.5,.1],color=:black,linewidth=1,label=false)
plot!([.5,0],[.5,.5],color=:black,linewidth=1,ls =:dash, label=false)
plot!([.5,.5],[.5,0],color=:black,linewidth=1,ls =:dash, label=false)
annotate!((-0.03,.47), "xᵢ",  9)
annotate!((.47,-0.03), "xᵢ",  9)
plot!([.5,.65],[.5,.1],color=:black,linewidth=1,label=false, c=:black, ls = :dot)

annotate!((.6,.3), "a.",  9)
annotate!((.53,.3), "b.",  9)
png(p, "C:/Users/jeffrey_andrews/OneDrive/Documents/phd/output/fehr.png")