/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_putstr.c                                        :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/04 11:57:21 by igomez            #+#    #+#             */
/*   Updated: 2014/11/06 16:50:20 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "libft.h"
#include <unistd.h>

void	ft_putstr(const char *str)
{
	size_t	i;

	if (!str)
		write(1, "(null)", 6);
	else
	{
		i = 0;
		while (str[i])
			i++;
		write(1, str, i);
	}
}
