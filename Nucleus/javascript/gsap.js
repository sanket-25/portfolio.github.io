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
var w = document.documentElement.clientWidth || window.innerWidth;
if (w >= 300) {
    gsap.from(".right-about", {
        scrollTrigger: {
            trigger: ".left-about",
            // markers: true,
            scroller: "body",
            endTrigger: ".right-about", // Custom trigger set to .right-about
            end: "bottom bottom",
			// scrub: 5,
            pin: true
        }
    });

    gsap.from(".left-skills", {
        scrollTrigger: {
            trigger: ".right-skills",
            // markers: true,
            scroller: "body",
            start: "top 0%",
			endTrigger: ".left-skills", // Custom trigger set to .right-about
            end: "bottom bottom",
            // scrub: 5,
            pin: true
        }
    });

    gsap.from(".right-work", {
        scrollTrigger: {
            trigger: ".left-work",
            // markers: true,
            scroller: "body",
            start: "top 0%",
            endTrigger: ".right-work", // Custom trigger set to .right-about
            end: "bottom bottom",
            // scrub: 5,
            pin: true
        }
    });
} else {
    // Code for smaller screens
}



