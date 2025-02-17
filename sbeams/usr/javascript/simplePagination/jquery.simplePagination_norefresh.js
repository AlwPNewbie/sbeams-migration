(function($){

    var methods = {
        init: function(options) {
            var o = $.extend({
                items: 1,
                itemsOnPage: 1,
                rs_set_name: '',
                table_id: '',
                pages: 0,
                displayedPages: 5,
                edges: 2,
                rs_page_number: 0,
                hrefTextPrefix: 'https://db.systemsbiology.net/devZS/sbeams/cgi/PeptideAtlas/ReadQueryResult.cgi?',
                hrefTextSuffix: '',
                prevText: 'Prev',
                nextText: 'Next',
                ellipseText: '&hellip;',
                ellipsePageSet: true,
                cssStyle: 'light-theme',
                listStyle: '',
                labelMap: [],
                selectOnClick: true,
                nextAtFront: false,
                invertPageOrder: false,
                useStartEdge: true,
                useEndEdge: true,
                onPageClick: function(pageNumber, event) {},
                onInit: function() {
                    // Callback triggered immediately after initialization
                }
            }, options || {});

            var self = this;
            o.pages = o.pages ? o.pages : Math.ceil(o.items / o.itemsOnPage) ? Math.ceil(o.items / o.itemsOnPage) : 1;
            if (o.rs_page_number)
                o.rs_page_number = o.rs_page_number - 1;
            else
                o.rs_page_number = !o.invertPageOrder ? 0 : o.pages - 1;
            o.halfDisplayed = o.displayedPages / 2;
            this.each(function() {
                self.addClass(o.cssStyle + ' simple-pagination').data('pagination', o);
                methods._draw.call(self);
            });
            o.onInit();
            return this;
        },

        selectPage: function(page) {
            methods._selectPage.call(this, page - 1);
            return this;
        },

        prevPage: function() {
            var o = this.data('pagination');
            if (!o.invertPageOrder) {
                if (o.rs_page_number > 0) {
                    methods._selectPage.call(this, o.rs_page_number - 1);
                }
            } else {
                if (o.rs_page_number < o.pages - 1) {
                    methods._selectPage.call(this, o.rs_page_number + 1);
                }
            }
            return this;
        },

        nextPage: function() {
            var o = this.data('pagination');
            if (!o.invertPageOrder) {
                if (o.rs_page_number < o.pages - 1) {
                    methods._selectPage.call(this, o.rs_page_number + 1);
                }
            } else {
                if (o.rs_page_number > 0) {
                    methods._selectPage.call(this, o.rs_page_number - 1);
                }
            }
            return this;
        },

        getPagesCount: function() {
            return this.data('pagination').pages;
        },

        setPagesCount: function(count) {
            this.data('pagination').pages = count;
        },

        getCurrentPage: function () {
            return this.data('pagination').rs_page_number + 1;
        },

        destroy: function(){
            this.empty();
            return this;
        },

        drawPage: function (page) {
            var o = this.data('pagination');
            o.rs_page_number = page - 1;
            this.data('pagination', o);
            methods._draw.call(this);
            return this;
        },

        redraw: function(){
            methods._draw.call(this);
            return this;
        },

        disable: function(){
            var o = this.data('pagination');
            o.disabled = true;
            this.data('pagination', o);
            methods._draw.call(this);
            return this;
        },

        enable: function(){
            var o = this.data('pagination');
            o.disabled = false;
            this.data('pagination', o);
            methods._draw.call(this);
            return this;
        },

        updateItems: function (newItems) {
            var o = this.data('pagination');
            o.items = newItems;
            o.pages = methods._getPages(o);
            this.data('pagination', o);
            methods._draw.call(this);
        },

        updateItemsOnPage: function (itemsOnPage) {
            var o = this.data('pagination');
            o.itemsOnPage = itemsOnPage;
            o.pages = methods._getPages(o);
            this.data('pagination', o);
            methods._selectPage.call(this, 0);
            return this;
        },

        getItemsOnPage: function() {
            return this.data('pagination').itemsOnPage;
        },

        _draw: function() {
            var	o = this.data('pagination'),
                interval = methods._getInterval(o),
                i,
                tagName;

            methods.destroy.call(this);

            tagName = (typeof this.prop === 'function') ? this.prop('tagName') : this.attr('tagName');

            var $panel = tagName === 'UL' ? this : $('<ul' + (o.listStyle ? ' class="' + o.listStyle + '"' : '') + '></ul>').appendTo(this);

            // Generate Prev link
            if (o.prevText) {
                methods._appendItem.call(this, !o.invertPageOrder ? o.rs_page_number - 1 : o.rs_page_number + 1, {text: o.prevText, classes: 'prev'});
            }

            // Generate Next link (if option set for at front)
            if (o.nextText && o.nextAtFront) {
                methods._appendItem.call(this, !o.invertPageOrder ? o.rs_page_number + 1 : o.rs_page_number - 1, {text: o.nextText, classes: 'next'});
            }

            // Generate start edges
            if (!o.invertPageOrder) {
                if (interval.start > 0 && o.edges > 0) {
                    if(o.useStartEdge) {
                        var end = Math.min(o.edges, interval.start);
                        for (i = 0; i < end; i++) {
                            methods._appendItem.call(this, i);
                        }
                    }
                    if (o.edges < interval.start && (interval.start - o.edges != 1)) {
                        $panel.append('<li class="disabled"><span class="ellipse">' + o.ellipseText + '</span></li>');
                    } else if (interval.start - o.edges == 1) {
                        methods._appendItem.call(this, o.edges);
                    }
                }
            } else {
                if (interval.end < o.pages && o.edges > 0) {
                    if(o.useStartEdge) {
                        var begin = Math.max(o.pages - o.edges, interval.end);
                        for (i = o.pages - 1; i >= begin; i--) {
                            methods._appendItem.call(this, i);
                        }
                    }

                    if (o.pages - o.edges > interval.end && (o.pages - o.edges - interval.end != 1)) {
                        $panel.append('<li class="disabled"><span class="ellipse">' + o.ellipseText + '</span></li>');
                    } else if (o.pages - o.edges - interval.end == 1) {
                        methods._appendItem.call(this, interval.end);
                    }
                }
            }

            // Generate interval links
            if (!o.invertPageOrder) {
                for (i = interval.start; i < interval.end; i++) {
                    methods._appendItem.call(this, i);
                }
            } else {
                for (i = interval.end - 1; i >= interval.start; i--) {
                    methods._appendItem.call(this, i);
                }
            }

            // Generate end edges
            if (!o.invertPageOrder) {
                if (interval.end < o.pages && o.edges > 0) {
                    if (o.pages - o.edges > interval.end && (o.pages - o.edges - interval.end != 1)) {
                        $panel.append('<li class="disabled"><span class="ellipse">' + o.ellipseText + '</span></li>');
                    } else if (o.pages - o.edges - interval.end == 1) {
                        methods._appendItem.call(this, interval.end);
                    }
                    if(o.useEndEdge) {
                        var begin = Math.max(o.pages - o.edges, interval.end);
                        for (i = begin; i < o.pages; i++) {
                            methods._appendItem.call(this, i);
                        }
                    }
                }
            } else {
                if (interval.start > 0 && o.edges > 0) {
                    if (o.edges < interval.start && (interval.start - o.edges != 1)) {
                        $panel.append('<li class="disabled"><span class="ellipse">' + o.ellipseText + '</span></li>');
                    } else if (interval.start - o.edges == 1) {
                        methods._appendItem.call(this, o.edges);
                    }

                    if(o.useEndEdge) {
                        var end = Math.min(o.edges, interval.start);
                        for (i = end - 1; i >= 0; i--) {
                            methods._appendItem.call(this, i);
                        }
                    }
                }
            }

            // Generate Next link (unless option is set for at front)
            if (o.nextText && !o.nextAtFront) {
                methods._appendItem.call(this, !o.invertPageOrder ? o.rs_page_number + 1 : o.rs_page_number - 1, {text: o.nextText, classes: 'next'});
            }
            if (o.ellipsePageSet && !o.disabled) {
                methods._ellipseClick.call(this, $panel);
            }

        },

        _getPages: function(o) {
            var pages = Math.ceil(o.items / o.itemsOnPage);
            return pages || 1;
        },

        _getInterval: function(o) {
            return {
                start: Math.ceil(o.rs_page_number > o.halfDisplayed ? Math.max(Math.min(o.rs_page_number - o.halfDisplayed, (o.pages - o.displayedPages)), 0) : 0),
                end: Math.ceil(o.rs_page_number > o.halfDisplayed ? Math.min(o.rs_page_number + o.halfDisplayed, o.pages) : Math.min(o.displayedPages, o.pages))
            };
        },

        _appendItem: function(pageIndex, opts) {
            var self = this, options, $link, o = self.data('pagination'), $linkWrapper = $('<li></li>'), $ul = self.find('ul');

            pageIndex = pageIndex < 0 ? 0 : (pageIndex < o.pages ? pageIndex : o.pages - 1);

            options = {
                text: pageIndex + 1,
                classes: ''
            };

            if (o.labelMap.length && o.labelMap[pageIndex]) {
                options.text = o.labelMap[pageIndex];
            }

            options = $.extend(options, opts || {});

            if (pageIndex == o.rs_page_number || o.disabled) {
                if (o.disabled || options.classes === 'prev' || options.classes === 'next') {
                    $linkWrapper.addClass('disabled');
                } else {
                    $linkWrapper.addClass('active');
                }
                $link = $('<span class="current">' + (options.text) + '</span>');

            } else {
                var table_str = ''
                if (o.table_id !== ''){
                    table_str = '&TABLE_id=' + o.table_id;
                }

                $link = $('<a href="javascript:void(0);" class="page-link">' + (options.text) + '</a>');
                $link.click(function(event){

                var href = o.hrefTextPrefix + 'rs_page=' + (pageIndex + 1) + '&rs_page_size=' + o.itemsOnPage +
                        '&rs_set_name=' + o.rs_set_name + table_str + o.hrefTextSuffix;

                    methods._selectPage.call(self, pageIndex, event, href);
                });
            }

            if (options.classes) {
                $link.addClass(options.classes);
            }

            $linkWrapper.append($link);

            if ($ul.length) {
                $ul.append($linkWrapper);
            } else {
                self.append($linkWrapper);
            }
        },

        _selectPage: function(pageIndex, event, href) {
            // Ensure that 'this' refers to a jQuery object
            var $this = this instanceof jQuery ? this : $(this);
            var o = $this.data('pagination');
						// Check if pagination options are available
						if (!o) {
								alert("Pagination plugin not initialized on this element.");
								return;
						}
            o.rs_page_number = pageIndex;
            if (o.selectOnClick) {
                $.ajax({
                    url: href,
                    method: 'GET',
                    dataType: 'json',
                    success: function(data) {
                        // Update the table with the fetched data
												updateTableWithJson(data, o.table_id);
                        $this.pagination('redraw', pageIndex);
                        //methods._draw.call(this);
                    },
                    error: function() {
                        console.error('Failed to fetch data from', href);
                    }
                });
            }
            return o.onPageClick(pageIndex + 1, event);
        },

        _ellipseClick: function($panel) {
            var self = this,
                o = this.data('pagination'),
                $ellip = $panel.find('.ellipse');
            $ellip.addClass('clickable').parent().removeClass('disabled');
            $ellip.click(function(event) {
                if (!o.disable) {
                    var $this = $(this),
                        val = (parseInt($this.parent().prev().text(), 10) || 0) + 1;
                    $this
                        .html('<input type="number" min="1" max="' + o.pages + '" step="1" value="' + val + '">')
                        .find('input')
                        .focus()
                        .click(function(event) {
                            // prevent input number arrows from bubbling a click event on $ellip
                            event.stopPropagation();
                        })
                        .keypress(function(event) {
                            var val = $(this).val();
                            if (event.which === 13 && val !== '') {
                                // enter to accept
                                if ((val > 0) && (val <= o.pages)){
                                    var href = o.hrefTextPrefix + 'rs_page=' + val + '&rs_page_size=' + o.itemsOnPage +
                                        '&rs_set_name=' + o.rs_set_name;
                                    if (o.table_id !== '') {
                                        href += '&TABLE_ID=' + o.table_id;
                                    }
                                    methods._selectPage.call(self, val - 1, event, href);
                                }
                            } else if (event.which === 27) {
                                // escape to cancel
                                $ellip.empty().html(o.ellipseText);
                            }
                        })
                        .bind('blur', function(event) {
                            var val = $(this).val();
                            if (val !== '') {
                                methods._selectPage.call(self, val - 1);
                            }
                            $ellip.empty().html(o.ellipseText);
                            return false;
                        });
                }
                return false;
            });
        }

    };

    $.fn.pagination = function(method) {
        // Method calling logic
        if (methods[method] && method.charAt(0) != '_') {
            return methods[method].apply(this, Array.prototype.slice.call(arguments, 1));
        } else if (typeof method === 'object' || ! method) {
            return methods.init.apply(this, arguments);
        } else {
            $.error('Method ' +  method + ' does not exist on jQuery.pagination');
        }

    };
})(jQuery);
function updateTableWithJson(jsonData, tableId) {
    var table = document.getElementById(tableId);
		var rows = table.querySelectorAll('tr');
		var rowsArray = Array.from(rows).slice(1);  // Exclude the header row


    // Define the background colors
    var colors = ["#f3f1e4", "#d3d1c4"];

    if (table){
			// Loop through each entry in jsonData
			jsonData.forEach(function(rowData, row_idx) {
				// Determine the background color based on the row index
				var bgColor = colors[Math.floor(row_idx / 3) % colors.length];

				// Check if a row exists for the current index
				if (row_idx < rowsArray.length) {
					// Update the existing row
					var row = rowsArray[row_idx];
					var cells = row.querySelectorAll('td');
					cells.forEach(function(cell, col_idx) {
            if (containsHttpLink(rowData[col_idx])){
              cell.innerHTML =  rowData[col_idx];
            }else{
						  cell.textContent = rowData[col_idx];
            }
					  cell.style.backgroundColor = bgColor; // Set the background color
					});
				} else {
					// Clone the last row if it doesn't exist
					var lastRow = rowsArray[rowsArray.length - 1];
					var newRow = lastRow.cloneNode(true);

					// Insert the cloned row into the table
					table.appendChild(newRow);
					var newCells = newRow.querySelectorAll('td');
					newCells.forEach(function(cell, col_idx) {
            if (containsHttpLink(rowData[col_idx])){
              cell.innerHTML =  rowData[col_idx];
            }else{
						  cell.textContent = rowData[col_idx];
            }
					  cell.style.backgroundColor = bgColor; // Set the background color
					});
				}
			});

			// Remove any extra rows if the original rows are more than jsonData length
			while (rowsArray.length > jsonData.length) {
				table.deleteRow(rowsArray.length);
				rowsArray.pop();
			}
		}
}
// Function to check if a string contains a URL starting with http:// or https://
function containsHttpLink(data) {
    // Regular expression to match http:// or https:// links
    const urlRegex = /(http:\/\/|https:\/\/)\S+/i;
    return urlRegex.test(data);
}
