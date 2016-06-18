/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_tabdel.c                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/01/09 16:44:38 by igomez            #+#    #+#             */
/*   Updated: 2015/01/09 16:44:39 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"
#include <stdlib.h>

void	ft_tabdel(char ***tab_ptr)
{
	size_t	i;

	if (!tab_ptr || !*tab_ptr)
		return ;
	i = 0;
	while ((*tab_ptr)[i])
		ft_strdel(*tab_ptr + i++);
	free(*tab_ptr);
	*tab_ptr = NULL;
}
