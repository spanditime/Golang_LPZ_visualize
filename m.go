package main

import (
	"fmt"
	"image/color"
	"math"
	"sort"

	"github.com/fogleman/gg"
)

type point struct {
	x float64
	y float64
}

func Newpoint(x, y float64) point {
	p := point{}
	p.x = x
	p.y = y
	return p
}

type line struct {
	ab point
	c  float64
}

func Newline(p1, p2 point) line {
	var ret line
	if p2.x == p1.x {
		ret.ab.x, ret.ab.y, ret.c = 1, 0, p1.x
	} else if p2.y == p1.y {
		ret.ab.x, ret.ab.y, ret.c = 0, 1, p1.y
	} else {
		ret.ab.x = (p2.x - p1.x)
		ret.ab.y = (p1.y - p2.y)
		ret.c = p1.y/ret.ab.y + p1.x/ret.ab.x
		ret.ab.x = 1 / ret.ab.x
		ret.ab.y = 1 / ret.ab.y
	}
	return ret
}

func IsNum(v float64) bool {
	return !math.IsNaN(v) && !math.IsInf(v, 0)
}

func (p *point) Add(v point) {
	p.x += v.x
	p.y += v.y
}

func (p *point) Div(v float64) {
	p.x /= v
	p.y /= v
}

func (p *point) Mul(v float64) {
	p.x *= v
	p.y *= v
}

func (p *point) Sub(v point) {
	p.x -= v.x
	p.y -= v.y
}

func (p point) Len() float64 {
	return math.Sqrt(p.x*p.x + p.y*p.y)
}

func DotProd(p1, p2 point) float64 {
	return p1.x*p2.x + p1.y*p2.y
}

// project a onto b
func Project(b, a point) point {
	b.Mul(DotProd(a, b) / DotProd(b, b))
	return b
}

func intersection(l1, l2 line) point {
	del := l1.ab.x*l2.ab.y - l2.ab.x*l1.ab.y
	xc := (l2.ab.y*l1.c - l1.ab.y*l2.c)
	yc := (l1.ab.x*l2.c - l2.ab.x*l1.c)
	if del != 0 && IsNum(del) && IsNum(xc) && IsNum(yc) {
		return Newpoint(
			xc/del,
			yc/del,
		)
	} else if IsNum(xc) && IsNum(yc) {
		return point{xc, yc}
	} else {
		return point{0, 0}
	}
}

type L struct {
	l  line
	op bool
}

func (l line) calculate(p point) float64 {
	return l.ab.x*p.x + l.ab.y*p.y
}

func (l L) calculate(p point) bool {
	if l.op {
		return l.l.calculate(p) >= l.l.c
	}
	return l.l.calculate(p) <= l.l.c
}

type LPZ struct {
	A   []L
	ODZ []L
}

type sorter struct {
	center point
}

func getAngle(center, p point) float64 {
	p.Sub(center)
	p.Div(p.Len())
	center.x = 1
	center.y = 0
	return math.Acos(DotProd(center, p) / (p.Len() * center.Len()))
}

func (lpz LPZ) getLines(mx, my float64) []struct{ p1, p2 point } {
	var ret []struct{ p1, p2 point }
	var lines []line = []line{
		{Newpoint(1, 0), 0},
		{Newpoint(1, 0), mx},
		{Newpoint(0, 1), 0},
		{Newpoint(0, 1), my},
	}
	for _, L := range lpz.A {
		var intersects []point
		for _, o := range lines {
			intersects = append(intersects, intersection(L.l, o))
		}
		i := 0
		for _, inter := range intersects {
			if inter.x >= 0 && inter.x <= mx && inter.y >= 0 && inter.y <= my {
				intersects[i] = inter
				i += 1
				if i == 2 {
					break
				}
			}
		}
		ret = append(ret, struct {
			p1 point
			p2 point
		}{p1: intersects[0], p2: intersects[1]})
	}
	return ret
}

func (lpz LPZ) getArea() ([]struct{ p1, p2 point }, []point, float64, float64) {
	points := []point{}
	srt := sorter{point{0, 0}}
	for i := range lpz.A {
		c := intersection(lpz.A[i].l, lpz.A[(i+1)%len(lpz.A)].l)
		points = append(points, c)
		srt.center.Add(c)
	}
	srt.center.Div(float64(len(points)))
	//arrange points
	sort.Slice(points, func(i, j int) bool {
		return getAngle(srt.center, points[i]) < getAngle(srt.center, points[j])
	})
	for _, l := range lpz.ODZ {
		var beg, end int = -1, -1
		for i := range points {
			var p point = points[i]
			cond := l.calculate(p)
			if beg == -1 && !cond {
				beg = i
			} else if beg != -1 && end == -1 && cond {
				end = i
			} else if beg != -1 && end != -1 && !cond {
				end = beg
				beg = i
				break
			}
		}
		if beg == -1 {
			fmt.Println("Error")
			continue
		}
		if end == -1 {
			end = 0
		}

		var bp1, bp2, ep1, ep2 point
		bp2, ep2 = points[beg], points[end]
		if beg == 0 {
			bp1 = points[len(points)-1]
		} else {
			bp1 = points[beg-1]
		}

		if end == 0 {
			ep1 = points[len(points)-1]
		} else {
			ep1 = points[end-1]
		}
		Line := l.l
		var p1 point = intersection(Line, Newline(bp1, bp2))
		var p2 point = intersection(Line, Newline(ep1, ep2))
		if end < beg {
			points = append(append([]point{p2}, points[end:beg]...), p1)
		} else {
			points = append(append(points[:beg], p1, p2), points[end:]...)
		}
		var i int = 0
		for i != len(points) {
			if points[i] == points[(i+1)%len(points)] {
				points = append(points[:i], points[i+1:]...)
			} else {
				i += 1
			}
		}
	}
	var my, mx float64 = 0., 0.
	for _, p := range points {
		my = math.Max(p.y, my)
		mx = math.Max(p.x, mx)
	}
	my, mx = math.Floor(my+3), math.Floor(mx+3)
	return lpz.getLines(mx, my), points, my, mx
}

type context struct {
	width  float64
	height float64
	ratio  float64
	off    float64
	xgz    bool
	ygz    bool
}

func (c context) Transform(x, y float64) (float64, float64) {
	var xr, yr float64
	if c.xgz {
		xr = c.off + x*c.ratio
	} else {
		xr = c.width - c.off - x*c.ratio
	}
	if c.ygz {
		yr = c.height - c.off - y*c.ratio
	} else {
		yr = c.off + y*c.ratio
	}
	return xr, yr
}

func (c context) TransformLine(x0, y0, x1, y1 float64) (float64, float64, float64, float64) {
	x0, y0 = c.Transform(x0, y0)
	x1, y1 = c.Transform(x1, y1)
	return x0, y0, x1, y1
}

// Z = 3x1 + x2
// -x1 + x2 >= 1
// x1 + 3x2 <= 15
// -2x1 + x2 <= 4
func main() {
	var lpz LPZ = LPZ{
		[]L{
			{line{Newpoint(-1, 1), 1}, true},
			{line{Newpoint(1, 3), 15}, false},
			{line{Newpoint(-2, 1), 4}, false},
		},
		[]L{
			{line{Newpoint(1, 0), 0}, true},
			{line{Newpoint(0, 1), 0}, true},
		},
	}
	var Z line = line{Newpoint(3, 1), 0}

	// var lpz LPZ = LPZ{
	// 	[]L{
	// 		{line{Newpoint(3, 2), 12}, false},
	// 		{line{Newpoint(1, 2), 4}, true},
	// 		{line{Newpoint(2, -1), 1}, false},
	// 	},
	// 	[]L{
	// 		{line{Newpoint(1, 0), 0}, true},
	// 		{line{Newpoint(0, 1), 0}, true},
	// 	},
	// }
	// var Z line = line{Newpoint(5, 1), 0}

	LINES, POINTS, mx, my := lpz.getArea()

	var MaxFP, MinFP point = POINTS[0], POINTS[0]
	var MaxF, MinF float64 = Z.calculate(MaxFP), Z.calculate(MinFP)
	for _, p := range POINTS {
		val := Z.calculate(p)
		if val > MaxF {
			MaxFP = p
			MaxF = val
		} else if val < MinF {
			MinFP = p
			MinF = val
		}
	}
	var GradMin, GradMax point = Project(Z.ab, MinFP), Project(Z.ab, MaxFP)

	mx = math.Max(mx, my)
	my = mx

	var ctx context = context{}

	ctx.width = 1000
	ctx.height = 1000
	ctx.off = 50
	ctx.ratio = (ctx.width - 2*ctx.off) / math.Max(mx, my)
	ctx.xgz = lpz.ODZ[0].op
	ctx.ygz = lpz.ODZ[1].op
	dc := gg.NewContext(int(ctx.width), int(ctx.height))

	//draw
	dc.SetStrokeStyle(gg.NewSolidPattern(color.White))
	dc.SetLineWidth(2)
	dc.DrawLine(ctx.TransformLine(0, 0, mx, 0))
	dc.DrawLine(ctx.TransformLine(0, 0, 0, my))
	dc.Stroke()
	dc.SetLineWidth(1)
	for i := 1; i <= int(math.Max(mx, my)); i++ {
		dc.DrawLine(ctx.TransformLine(0, float64(i), mx, float64(i)))
		dc.DrawLine(ctx.TransformLine(float64(i), 0, float64(i), my))
	}
	dc.Stroke()
	//draw area
	grad := gg.NewLinearGradient(
		ctx.off+GradMin.x*ctx.ratio,
		ctx.height-ctx.off-GradMin.y*ctx.ratio,
		ctx.off+GradMax.x*ctx.ratio,
		ctx.height-ctx.off-GradMax.y*ctx.ratio,
	)
	grad.AddColorStop(0, color.RGBA{220, 8, 7, 255})
	grad.AddColorStop(1, color.RGBA{10, 232, 11, 255})
	dc.SetStrokeStyle(grad)
	dc.SetFillStyle(grad)
	dc.SetLineWidth(4)
	for i := range POINTS {
		var p point = POINTS[(i+1)%len(POINTS)]
		if i == 0 {
			dc.MoveTo(ctx.Transform(POINTS[0].x, POINTS[0].y))
		}
		dc.LineTo(ctx.Transform(p.x, p.y))
	}
	dc.Fill()
	dc.SetStrokeStyle(grad)
	dc.SetLineWidth(4)
	dc.MoveTo(ctx.Transform(GradMin.x, GradMin.y))
	dc.LineTo(ctx.Transform(GradMax.x, GradMax.y))
	dc.Stroke()
	//drawLines

	dc.SetStrokeStyle(gg.NewSolidPattern(color.White))
	dc.SetLineWidth(2)
	for _, line := range LINES {
		dc.DrawLine(
			ctx.TransformLine(line.p1.x, line.p1.y, line.p2.x, line.p2.y),
		)
		dc.Stroke()
	}
	dc.SetLineWidth(1)
	dc.DrawLine(ctx.TransformLine(GradMax.x, GradMax.y, MaxFP.x, MaxFP.y))
	dc.DrawLine(ctx.TransformLine(GradMin.x, GradMin.y, MinFP.x, MinFP.y))
	dc.Stroke()
	{
		x, y := ctx.Transform(MinFP.x, MinFP.y)
		dc.SetFillStyle(gg.NewSolidPattern(color.RGBA{220, 8, 7, 255}))
		dc.DrawCircle(x, y, 5)
		dc.Fill()

		x, y = ctx.Transform(MaxFP.x, MaxFP.y)
		dc.SetFillStyle(gg.NewSolidPattern(color.RGBA{10, 232, 11, 255}))
		dc.DrawCircle(x, y, 5)
	}
	dc.Fill()

	//save png
	dc.SavePNG("result.png")
}
