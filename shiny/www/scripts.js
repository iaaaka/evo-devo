// Start only when document is fully loaded
$(document).ready(function() {
	
	//$('#gene_title').hide();

	$('#click').click(function() {
		if ($('#gene').val() !== null ) {
			$('#update').trigger('click');
			$('#species').trigger('click');
			$('#gene').val(null);
			$('#gene_title').css('display', 'initial');
			return;
		}
	});

	setInterval(function(){
		if ($('html').attr('class')=='shiny-busy') {
			$('#to_hide').css('opacity', '0');
			$('div.pbusy').show();
			//console.log("busy");
		} else {
			$('div.pbusy').hide();
			$('#to_hide').delay(800).css('opacity', '1');
			//console.log("done");
		}
	}, 100);

	$(window).on('load resize', function() {
		var h = $('#absNav').height();
		$('.sideChoices').css('margin-top', h + 15);
		$('#link_lab').css('top', (h/2) - 10);
	});
	
	if ($('#gene_title')) {
	  $('#mainPicture').css('left', '28px');
	}

/*	$(document).on("keypress", function (e) {
		if (e.keyCode == 13) {
			if ($('#gene').val() !== null ) {
				$('#click').trigger('click');
			 // $('#gene').val(null);
				return;
			}
		}
	});
*/
});