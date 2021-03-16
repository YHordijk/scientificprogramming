$(document).ready(function() {
    // make citation labels clickable
    $('.docutils.citation').each(function() {
        var label = $(this).find('td.label');
        var temp = label.text();
        var a = $('<a/>', {
            text: temp,
            href: '#' + this.id
        });
        label.html(a);
    });
});