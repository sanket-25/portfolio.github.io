// gsap.from(".box", {
// 	scale:0,
// 	delay:1,
// 	duration:3,
// 	rotate:360,
// 	x: 990,
// 	scrollTrigger: {
// 		trigger: ".box",
// 		markers: true,
// 		scroller: "body",
// 		start: "top 50%",
// 		end: "top 30%",
// 		scrub: 5

// 	}
// })



gsap.from(".right-about", {
	scrollTrigger: {
		trigger: ".left-about",
		// markers: true,
		scroller: "body",
		start: "top 0%",
		end: "bottom 0%",
		// scrub: 5,
		pin: true

	}
})

gsap.from(".left-skills", {
	scrollTrigger: {
		trigger: ".right-skills",
		// markers: true,
		scroller: "body",
		start: "top 0%",
		end: "bottom 10%",
		// scrub: 5,
		pin: true

	}
})

gsap.from(".right-work", {
	scrollTrigger: {
		trigger: ".left-work",
		markers: true,
		scroller: "body",
		start: "top 0%",
		end: "bottom -180%",
		// scrub: 5,
		pin: true

	}
})