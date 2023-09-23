// $(document).ready(function() {
//     $("#color_mode").on("change", function () {
//         colorModePreview(this);
//     })
// });

// function colorModePreview(ele) {
//     if($(ele).prop("checked") == true){
//         $('body').addClass('dark-preview');
//         $('body').removeClass('white-preview');
//     }
//     else if($(ele).prop("checked") == false){
//         $('body').addClass('white-preview');
//         $('body').removeClass('dark-preview');
//     }
// }


gsap.from(".box", {
	scale:0,
	delay:1,
	duration:3,
	rotate:360,
	x: 990,
	scrollTrigger: {
		trigger: ".box",
		markers: true,
		scroller: "body",
		start: "top 50%",
		end: "top 30%",
		scrub: 5

	}
})