/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_putaddr.c                                       :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2014/11/18 09:53:54 by igomez            #+#    #+#             */
/*   Updated: 2014/11/28 12:40:49 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include <unistd.h>
#include "libft.h"

void	ft_putaddr(void *str)
{
	char				temp;
	char				*base_str;
	int					i;
	unsigned long int	n;
	unsigned long int	div;

	n = (unsigned long int)str;
	base_str = "0123456789abcdef";
	write(1, "0x", 2);
	i = 7;
	div = 1;
	while (!(n / div < 16))
	{
		div *= 16;
		i--;
	}
	while (i-- > 0)
		ft_putchar('0');
	while (div > 0)
	{
		temp = n / div;
		write(1, base_str + temp, 1);
		n %= div;
		div /= 16;
	}
}
